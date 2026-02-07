"""Joint strain-function EM inference engine for StrainSift.

THE CORE MODULE. Simultaneously estimates:
- pi: strain abundance vector (Dirichlet prior)
- phi: strain-function attribution matrix (Beta prior)
- gamma: per-read strain assignment posteriors

With mobile element pool, Fisher information diagnostics,
BIC model selection, multiple restarts, and ablation mode.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.special import logsumexp

from strainsift.types import (
    ConvergenceDiagnostics,
    ReadAssignment,
    ReadClassification,
    StrainProfile,
)

logger = logging.getLogger(__name__)


@dataclass
class EMConfig:
    """Configuration for EM inference.

    Attributes:
        max_iterations: Maximum EM iterations per restart.
        convergence_threshold: Relative change in log-likelihood for convergence.
        n_restarts: Number of random restarts.
        alpha: Dirichlet prior concentration (< 1 for sparsity).
        beta_a: Beta prior shape a for phi.
        beta_b: Beta prior shape b for phi.
        min_containment: Minimum containment to consider a candidate.
        mobile_pool: Whether to include the mobile element pool.
        fix_phi: If True, don't update phi (strain-only ablation mode).
        fix_pi: If True, don't update pi (function-only ablation mode).
        random_seed: Random seed for reproducibility.
    """

    max_iterations: int = 200
    convergence_threshold: float = 1e-6
    n_restarts: int = 10
    alpha: float = 0.5
    beta_a: float = 1.5
    beta_b: float = 1.5
    min_containment: float = 0.05
    mobile_pool: bool = True
    fix_phi: bool = False
    fix_pi: bool = False
    random_seed: int = 42


# Internal read representation for efficient EM computation
@dataclass
class _ReadData:
    """Preprocessed read data for EM.

    Attributes:
        strain_indices: Indices into strain_ids for candidate strains.
        gene_indices: Indices into gene_family_ids (-1 for intergenic/None).
        containments: Containment scores for each candidate.
    """

    strain_indices: list[int]
    gene_indices: list[int]
    containments: list[float]


class JointEM:
    """Joint strain-function EM inference engine.

    The core algorithm of StrainSift. Simultaneously estimates strain
    abundances (pi), gene-strain attribution (phi), and per-read strain
    assignments (gamma) via Expectation-Maximization with MAP estimation.

    The key insight: if a read hits a gene family g, strains known to
    carry g (high phi[k][g]) get higher responsibility. This is where
    "function helps resolve strains."

    Usage:
        em = JointEM(config)
        profile = em.fit(classifications, strain_ids, gene_family_ids)
    """

    def __init__(self, config: EMConfig | None = None):
        self.config = config or EMConfig()

    def fit(
        self,
        read_classifications: list[ReadClassification],
        strain_ids: list[str],
        gene_family_ids: list[str],
        reference_phi: np.ndarray | None = None,
    ) -> StrainProfile:
        """Run joint EM inference.

        Args:
            read_classifications: Per-read candidate lists from classification.
            strain_ids: List of K strain identifiers.
            gene_family_ids: List of G gene family identifiers.
            reference_phi: Optional K x G matrix of reference annotations
                (1.0 = strain carries gene, 0.0 = does not).

        Returns:
            StrainProfile with all estimates and uncertainty.
        """
        cfg = self.config
        K = len(strain_ids)
        G = len(gene_family_ids)

        # Preprocess reads into efficient internal format
        read_data = self._preprocess_reads(
            read_classifications, strain_ids, gene_family_ids
        )

        # Filter to classified reads only
        read_data = [rd for rd in read_data if rd.strain_indices]
        N = len(read_data)

        if N == 0:
            logger.warning("No classified reads for EM inference")
            return self._empty_profile(
                strain_ids, gene_family_ids, len(read_classifications)
            )

        logger.info(
            "EM inference: %d classified reads, %d strains, %d gene families",
            N, K, G,
        )

        # Number of components: K strains + optional mobile pool
        K_total = K + 1 if cfg.mobile_pool else K

        # Build Beta prior hyperparameter matrices from reference
        prior_a, prior_b = self._build_phi_priors(K, G, reference_phi)

        # Multiple restarts
        best_ll = -np.inf
        best_pi = None
        best_phi = None
        best_gamma = None
        best_ll_trace = []
        restart_lls = []

        base_rng = np.random.default_rng(cfg.random_seed)

        for restart in range(cfg.n_restarts):
            seed = base_rng.integers(0, 2**31)
            rng = np.random.default_rng(seed)

            # Initialize
            pi, phi = self._initialize(K, G, K_total, rng, reference_phi)

            # EM iterations
            ll_trace = []
            converged = False

            for iteration in range(cfg.max_iterations):
                # E-step
                gamma, ll = self._e_step(pi, phi, read_data, K, G, K_total)
                ll_trace.append(ll)

                # Check convergence
                if len(ll_trace) >= 2:
                    prev_ll = ll_trace[-2]
                    if abs(prev_ll) > 0:
                        rel_change = abs(ll - prev_ll) / abs(prev_ll)
                    else:
                        rel_change = abs(ll - prev_ll)
                    if rel_change < cfg.convergence_threshold:
                        converged = True
                        break

                # M-step
                pi, phi = self._m_step(
                    gamma, read_data, K, G, K_total, pi, phi, prior_a, prior_b
                )

            restart_lls.append(ll_trace[-1] if ll_trace else -np.inf)

            if ll_trace and ll_trace[-1] > best_ll:
                best_ll = ll_trace[-1]
                best_pi = pi.copy()
                best_phi = phi.copy()
                best_gamma = gamma.copy()
                best_ll_trace = list(ll_trace)

            logger.debug(
                "Restart %d/%d: LL=%.4f, %d iterations, converged=%s",
                restart + 1, cfg.n_restarts, ll_trace[-1] if ll_trace else -np.inf,
                len(ll_trace), converged,
            )

        # Use best restart results
        pi = best_pi
        phi = best_phi
        gamma = best_gamma

        # Normalize phi per strain for interpretable gene presence output.
        # During EM, phi represents a read-fraction rate (~ gene_length /
        # genome_length per present gene). Rescaling by the median of
        # positive values converts the rate into a gene presence indicator
        # that thresholds correctly at 0.5. Median is more robust than max
        # to uneven gene coverage.
        for k in range(K):
            positive = phi[k][phi[k] > 1e-4]
            if len(positive) > 0:
                norm_val = np.median(positive)
                if norm_val > 1e-6:
                    phi[k] = phi[k] / norm_val
        phi = np.clip(phi, 1e-6, 1.0 - 1e-6)

        # Convergence diagnostics
        multimodal = False
        if len(restart_lls) >= 2:
            ll_range = max(restart_lls) - min(restart_lls)
            if ll_range > 10:
                multimodal = True
                logger.warning(
                    "Multimodal EM: LL range across restarts = %.2f", ll_range
                )

        convergence = ConvergenceDiagnostics(
            log_likelihoods=best_ll_trace,
            n_restarts=cfg.n_restarts,
            restart_log_likelihoods=restart_lls,
            best_restart=int(np.argmax(restart_lls)),
            converged=len(best_ll_trace) < cfg.max_iterations,
            n_iterations=len(best_ll_trace),
            multimodal_warning=multimodal,
        )

        # Fisher information
        fisher_info = self._compute_fisher_information(
            pi, phi, gamma, read_data, K, K_total
        )

        # Confidence intervals
        pi_ci_lower, pi_ci_upper, phi_ci_lower, phi_ci_upper = (
            self._compute_confidence_intervals(
                pi, phi, fisher_info, gamma, read_data,
                K, G, K_total, N, prior_a, prior_b,
            )
        )

        # BIC
        bic = self._compute_bic(best_ll, K, G, N, cfg.mobile_pool)

        # Identifiability scores from Fisher information
        ident_scores = self._compute_identifiability_scores(
            fisher_info, strain_ids, K
        )

        # Build per-read assignments
        read_assignments = self._build_read_assignments(
            gamma, read_data, read_classifications, strain_ids,
            gene_family_ids, K, K_total,
        )

        # Mobile pool proportion
        mobile_prop = float(pi[K]) if cfg.mobile_pool and len(pi) > K else 0.0

        species_id = "unknown"
        # Try to infer species from strain IDs
        if strain_ids:
            parts = strain_ids[0].rsplit("_", 1)
            if len(parts) > 1:
                species_id = "inferred_species"

        return StrainProfile(
            species_id=species_id,
            strain_ids=list(strain_ids),
            pi=pi,
            pi_ci_lower=pi_ci_lower,
            pi_ci_upper=pi_ci_upper,
            phi=phi,
            phi_ci_lower=phi_ci_lower,
            phi_ci_upper=phi_ci_upper,
            gene_family_ids=list(gene_family_ids),
            mobile_pool_proportion=mobile_prop,
            fisher_information=fisher_info,
            identifiability_scores=ident_scores,
            convergence=convergence,
            bic=bic,
            n_reads_classified=N,
            n_reads_total=len(read_classifications),
            read_assignments=read_assignments,
        )

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def _initialize(
        self,
        K: int,
        G: int,
        K_total: int,
        rng: np.random.Generator,
        reference_phi: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Initialize pi and phi for one EM restart.

        Returns (pi, phi) where pi has length K_total and phi has shape (K, G).
        """
        cfg = self.config

        # Initialize pi from perturbed Dirichlet
        alpha_init = np.full(K_total, cfg.alpha + 0.5)
        pi = rng.dirichlet(alpha_init)

        # Initialize phi
        if reference_phi is not None:
            # Perturb reference phi slightly
            noise = rng.uniform(-0.1, 0.1, size=(K, G))
            phi = np.clip(reference_phi + noise, 0.01, 0.99)
        else:
            # Draw from Beta prior
            phi = rng.beta(cfg.beta_a, cfg.beta_b, size=(K, G))

        return pi, phi

    # ------------------------------------------------------------------
    # E-step
    # ------------------------------------------------------------------

    def _e_step(
        self,
        pi: np.ndarray,
        phi: np.ndarray,
        read_data: list[_ReadData],
        K: int,
        G: int,
        K_total: int,
    ) -> tuple[np.ndarray, float]:
        """E-step: compute responsibilities gamma and penalized log-likelihood.

        For each read r with candidates:
        log_gamma[r][k] = log(pi[k]) + log(C(read_r, strain_k)) + log(phi[k][f_r])

        Normalized via logsumexp. Mobile pool uses mean containment
        and uniform 1/G for phi.

        The returned log-likelihood includes prior terms (Dirichlet on pi,
        Beta on phi) so that the MAP-EM guarantees monotonic non-decrease.

        Returns:
            gamma: (N, K_total) responsibility matrix.
            penalized_ll: Penalized log-likelihood (data LL + log prior).
        """
        cfg = self.config
        N = len(read_data)
        gamma = np.zeros((N, K_total))
        data_ll = 0.0

        log_pi = np.log(np.maximum(pi, 1e-300))

        for r, rd in enumerate(read_data):
            if not rd.strain_indices:
                continue

            log_probs = np.full(K_total, -np.inf)

            for idx in range(len(rd.strain_indices)):
                k = rd.strain_indices[idx]
                c = rd.containments[idx]
                g = rd.gene_indices[idx]

                if k >= K:
                    continue

                log_c = np.log(max(c, 1e-300))

                if g >= 0 and g < G:
                    log_phi_kg = np.log(max(phi[k, g], 1e-300))
                else:
                    log_phi_kg = 0.0

                log_prob = log_pi[k] + log_c + log_phi_kg

                # Aggregate: take max if strain appears multiple times
                if log_prob > log_probs[k]:
                    log_probs[k] = log_prob

            # Mobile pool component
            if cfg.mobile_pool and K_total > K:
                mean_c = np.mean(rd.containments) if rd.containments else 0.01
                log_c_mobile = np.log(max(mean_c * 0.5, 1e-300))
                log_phi_mobile = -np.log(max(G, 1))
                log_probs[K] = log_pi[K] + log_c_mobile + log_phi_mobile

            # Normalize
            valid_mask = log_probs > -np.inf
            if not np.any(valid_mask):
                gamma[r, :] = 1.0 / K_total
                continue

            log_norm = logsumexp(log_probs[valid_mask])
            for k in range(K_total):
                if log_probs[k] > -np.inf:
                    gamma[r, k] = np.exp(log_probs[k] - log_norm)

            data_ll += log_norm

        # Add log prior terms for MAP-EM monotonicity guarantee
        # Dirichlet prior on pi
        log_prior = 0.0
        for k in range(K_total):
            if pi[k] > 1e-300:
                log_prior += (cfg.alpha - 1.0) * np.log(pi[k])

        # Beta prior on phi
        for k in range(K):
            for g in range(G):
                p = np.clip(phi[k, g], 1e-300, 1.0 - 1e-300)
                log_prior += (cfg.beta_a - 1.0) * np.log(p)
                log_prior += (cfg.beta_b - 1.0) * np.log(1.0 - p)

        penalized_ll = data_ll + log_prior
        return gamma, penalized_ll

    # ------------------------------------------------------------------
    # M-step
    # ------------------------------------------------------------------

    def _m_step(
        self,
        gamma: np.ndarray,
        read_data: list[_ReadData],
        K: int,
        G: int,
        K_total: int,
        pi: np.ndarray,
        phi: np.ndarray,
        prior_a: np.ndarray,
        prior_b: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """M-step: update pi and phi using MAP estimates.

        pi[k] = (sum_r gamma[r][k] + alpha - 1) / (N + K_total*(alpha-1))
        phi[k][g] = (sum_{r: f_r=g} gamma[r][k] + a_kg - 1) /
                    (sum_r gamma[r][k] + a_kg + b_kg - 2)

        Respects ablation controls.
        """
        cfg = self.config
        N = len(read_data)

        # Update pi (MAP Dirichlet)
        if not cfg.fix_pi:
            gamma_sum = gamma.sum(axis=0)  # (K_total,)
            new_pi = gamma_sum + cfg.alpha - 1.0
            new_pi = np.maximum(new_pi, 1e-10)  # avoid negative from prior
            new_pi /= new_pi.sum()
            pi = new_pi
        elif cfg.fix_pi:
            # Function-only mode: keep pi uniform
            pi = np.ones(K_total) / K_total

        # Update phi (MAP Beta)
        if not cfg.fix_phi:
            # Accumulate: for each strain k and gene g, sum gamma for reads
            # that map to gene g in strain k.
            # CRITICAL: the denominator uses only gamma from reads that hit
            # ANY gene for the strain (genic reads), not intergenic reads.
            # Otherwise intergenic reads massively dilute phi toward 0.
            gene_gamma = np.zeros((K, G))  # numerator: gamma for gene g
            strain_gamma_genic = np.zeros(K)  # denominator: gamma from genic reads

            for r, rd in enumerate(read_data):
                for idx in range(len(rd.strain_indices)):
                    k = rd.strain_indices[idx]
                    g = rd.gene_indices[idx]
                    if k < K and 0 <= g < G:
                        gene_gamma[k, g] += gamma[r, k]
                        strain_gamma_genic[k] += gamma[r, k]

            # MAP Beta update
            numerator = gene_gamma + prior_a - 1.0
            denominator = strain_gamma_genic[:, np.newaxis] + prior_a + prior_b - 2.0

            # Clip to avoid division by zero and ensure valid probabilities
            numerator = np.maximum(numerator, 1e-10)
            denominator = np.maximum(denominator, numerator + 1e-10)

            phi = numerator / denominator
            phi = np.clip(phi, 1e-6, 1.0 - 1e-6)

        return pi, phi

    # ------------------------------------------------------------------
    # Fisher information and diagnostics
    # ------------------------------------------------------------------

    def _compute_fisher_information(
        self,
        pi: np.ndarray,
        phi: np.ndarray,
        gamma: np.ndarray,
        read_data: list[_ReadData],
        K: int,
        K_total: int,
    ) -> np.ndarray:
        """Compute observed Fisher information matrix for strain proportions.

        F[j,k] = sum_r gamma[r][j] * (delta(j,k)/pi[j] - 1) / pi[j]

        Used for confidence intervals and identifiability diagnostics.
        Returns K x K Fisher information matrix (for strains only, not mobile pool).
        """
        N = len(read_data)
        fisher = np.zeros((K, K))

        if N == 0:
            return fisher

        # Observed Fisher information from sufficient statistics
        gamma_strains = gamma[:, :K]  # (N, K)
        n_k = gamma_strains.sum(axis=0)  # effective counts per strain

        for j in range(K):
            for k in range(K):
                if j == k:
                    # Diagonal: sum_r gamma[r][j] / pi[j]^2
                    if pi[j] > 1e-10:
                        fisher[j, k] = n_k[j] / (pi[j] ** 2)
                else:
                    # Off-diagonal: -sum_r gamma[r][j]*gamma[r][k] / (pi[j]*pi[k])
                    # Approximated by the multinomial Fisher information
                    if pi[j] > 1e-10 and pi[k] > 1e-10:
                        # Cross-term from multinomial
                        fisher[j, k] = 0.0  # multinomial Fisher is diagonal

        # Add regularization for numerical stability
        fisher += np.eye(K) * 1e-6

        return fisher

    def _compute_identifiability_scores(
        self,
        fisher_info: np.ndarray,
        strain_ids: list[str],
        K: int,
    ) -> dict[tuple[str, str], float]:
        """Compute pairwise identifiability scores from Fisher information.

        High score = strains are distinguishable. Low score = strains may
        need to be merged.
        """
        scores = {}
        if K < 2:
            return scores

        # Use eigenvalues of the Fisher information matrix
        try:
            eigvals = np.linalg.eigvalsh(fisher_info)
            min_eigval = float(np.min(eigvals))
        except np.linalg.LinAlgError:
            min_eigval = 0.0

        for i in range(K):
            for j in range(i + 1, K):
                # Pairwise score: based on difference in Fisher info diagonals
                # and the minimum eigenvalue
                score = min(fisher_info[i, i], fisher_info[j, j])
                score = max(score, min_eigval)
                scores[(strain_ids[i], strain_ids[j])] = float(score)

        return scores

    # ------------------------------------------------------------------
    # Confidence intervals
    # ------------------------------------------------------------------

    def _compute_confidence_intervals(
        self,
        pi: np.ndarray,
        phi: np.ndarray,
        fisher_info: np.ndarray,
        gamma: np.ndarray,
        read_data: list[_ReadData],
        K: int,
        G: int,
        K_total: int,
        N: int,
        prior_a: np.ndarray,
        prior_b: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute 95% CIs for pi and phi.

        For pi: from inverse Fisher information.
        For phi: from posterior Beta distribution.
        """
        z = 1.96  # 95% CI

        # Pi CIs from Fisher information
        pi_ci_lower = np.zeros(K_total)
        pi_ci_upper = np.ones(K_total)

        try:
            # Invert Fisher information for variance
            if K > 0 and np.linalg.det(fisher_info) > 1e-20:
                cov = np.linalg.inv(fisher_info)
                pi_var = np.diag(cov)
                pi_var = np.maximum(pi_var, 0.0)  # numerical safety
                pi_se = np.sqrt(pi_var)

                for k in range(K):
                    pi_ci_lower[k] = max(0.0, pi[k] - z * pi_se[k])
                    pi_ci_upper[k] = min(1.0, pi[k] + z * pi_se[k])

                # Mobile pool CI (simple heuristic)
                if K_total > K:
                    mobile_se = np.sqrt(pi[K] * (1 - pi[K]) / max(N, 1))
                    pi_ci_lower[K] = max(0.0, pi[K] - z * mobile_se)
                    pi_ci_upper[K] = min(1.0, pi[K] + z * mobile_se)
            else:
                # Fisher info not invertible — use wide CIs
                for k in range(K_total):
                    pi_se = np.sqrt(pi[k] * (1 - pi[k]) / max(N, 1))
                    pi_ci_lower[k] = max(0.0, pi[k] - z * pi_se)
                    pi_ci_upper[k] = min(1.0, pi[k] + z * pi_se)
        except np.linalg.LinAlgError:
            for k in range(K_total):
                pi_se = np.sqrt(pi[k] * (1 - pi[k]) / max(N, 1))
                pi_ci_lower[k] = max(0.0, pi[k] - z * pi_se)
                pi_ci_upper[k] = min(1.0, pi[k] + z * pi_se)

        # Phi CIs from posterior Beta
        phi_ci_lower = np.zeros((K, G))
        phi_ci_upper = np.ones((K, G))

        # Compute effective counts per strain-gene pair
        gene_gamma = np.zeros((K, G))
        strain_gamma = gamma[:, :K].sum(axis=0)

        for r, rd in enumerate(read_data):
            for idx in range(len(rd.strain_indices)):
                k = rd.strain_indices[idx]
                g = rd.gene_indices[idx]
                if k < K and 0 <= g < G:
                    gene_gamma[k, g] += gamma[r, k]

        # Posterior Beta parameters
        post_a = gene_gamma + prior_a
        post_b = (strain_gamma[:, np.newaxis] - gene_gamma) + prior_b

        # Beta distribution CI
        from scipy.stats import beta as beta_dist

        for k in range(K):
            for g in range(G):
                a = max(post_a[k, g], 0.01)
                b = max(post_b[k, g], 0.01)
                try:
                    phi_ci_lower[k, g] = beta_dist.ppf(0.025, a, b)
                    phi_ci_upper[k, g] = beta_dist.ppf(0.975, a, b)
                except (ValueError, FloatingPointError):
                    phi_ci_lower[k, g] = 0.0
                    phi_ci_upper[k, g] = 1.0

        return pi_ci_lower, pi_ci_upper, phi_ci_lower, phi_ci_upper

    # ------------------------------------------------------------------
    # Model selection
    # ------------------------------------------------------------------

    def _compute_bic(
        self, log_likelihood: float, K: int, G: int, N: int, mobile_pool: bool
    ) -> float:
        """Bayesian Information Criterion.

        BIC = -2 * LL + n_params * ln(N)
        n_params = (K-1) for pi + K*G for phi + (1 if mobile pool)
        """
        n_params = (K - 1) + K * G
        if mobile_pool:
            n_params += 1  # mobile pool proportion
        if N <= 0:
            return np.inf
        return -2.0 * log_likelihood + n_params * np.log(N)

    def select_k(
        self,
        read_classifications: list[ReadClassification],
        strain_ids_per_k: dict[int, list[str]],
        gene_family_ids: list[str],
        k_max: int = 10,
        reference_phi_per_k: dict[int, np.ndarray] | None = None,
    ) -> tuple[int, StrainProfile]:
        """Model selection: find optimal K by BIC.

        Runs EM for K = 1, 2, ..., k_max.
        Additionally applies Fisher information merge criterion.

        Args:
            read_classifications: Classified reads.
            strain_ids_per_k: Dict mapping K -> list of strain IDs for that K.
            gene_family_ids: Gene family identifiers.
            k_max: Maximum K to try.
            reference_phi_per_k: Optional dict mapping K -> reference phi matrix.

        Returns:
            Tuple of (optimal_K, best_StrainProfile).
        """
        best_bic = np.inf
        best_k = 1
        best_profile = None

        for k_val in range(1, k_max + 1):
            if k_val not in strain_ids_per_k:
                continue

            strain_ids = strain_ids_per_k[k_val]
            ref_phi = None
            if reference_phi_per_k and k_val in reference_phi_per_k:
                ref_phi = reference_phi_per_k[k_val]

            profile = self.fit(
                read_classifications, strain_ids, gene_family_ids, ref_phi
            )

            logger.info("K=%d: BIC=%.2f, LL=%.2f", k_val, profile.bic,
                        profile.convergence.log_likelihoods[-1] if profile.convergence.log_likelihoods else float('nan'))

            if profile.bic < best_bic:
                best_bic = profile.bic
                best_k = k_val
                best_profile = profile

        if best_profile is None:
            # Fallback: use K=1
            strain_ids = strain_ids_per_k.get(1, ["strain_0"])
            best_profile = self.fit(
                read_classifications, strain_ids, gene_family_ids
            )
            best_k = 1

        logger.info("Model selection: optimal K=%d (BIC=%.2f)", best_k, best_bic)
        return best_k, best_profile

    # ------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------

    def _preprocess_reads(
        self,
        read_classifications: list[ReadClassification],
        strain_ids: list[str],
        gene_family_ids: list[str],
    ) -> list[_ReadData]:
        """Convert ReadClassification objects to efficient internal format."""
        strain_idx = {sid: i for i, sid in enumerate(strain_ids)}
        gene_idx = {gid: i for i, gid in enumerate(gene_family_ids)}

        result = []
        for rc in read_classifications:
            rd = _ReadData(
                strain_indices=[],
                gene_indices=[],
                containments=[],
            )

            if rc.is_classified:
                for cand in rc.candidates:
                    k = strain_idx.get(cand.strain_id)
                    if k is None:
                        continue
                    if cand.containment < self.config.min_containment:
                        continue

                    g = gene_idx.get(cand.gene_family_id, -1) if cand.gene_family_id else -1

                    rd.strain_indices.append(k)
                    rd.gene_indices.append(g)
                    rd.containments.append(cand.containment)

            result.append(rd)

        return result

    # ------------------------------------------------------------------
    # Output building
    # ------------------------------------------------------------------

    def _build_read_assignments(
        self,
        gamma: np.ndarray,
        read_data: list[_ReadData],
        read_classifications: list[ReadClassification],
        strain_ids: list[str],
        gene_family_ids: list[str],
        K: int,
        K_total: int,
    ) -> list[ReadAssignment]:
        """Build per-read ReadAssignment objects from gamma.

        gamma and read_data are aligned (same length, same ordering) —
        both correspond to the classified-and-filtered reads.  We walk
        through all read_classifications and match them to gamma rows
        by read_id.
        """
        # read_data and gamma are aligned: both are the result of
        # _preprocess_reads() filtered to entries with strain_indices.
        # We rebuild the mapping by replaying the same filtering on the
        # original read_classifications.

        # Step 1: recompute which reads survived preprocessing + filtering
        strain_idx_set = {sid for sid in strain_ids}
        gamma_map: dict[str, int] = {}
        gamma_idx = 0
        for rc in read_classifications:
            if not rc.is_classified:
                continue
            # Check if this read would have strain_indices after preprocessing
            has_valid = False
            for cand in rc.candidates:
                if (cand.strain_id in strain_idx_set
                        and cand.containment >= self.config.min_containment):
                    has_valid = True
                    break
            if has_valid and gamma_idx < len(read_data):
                gamma_map[rc.read_id] = gamma_idx
                gamma_idx += 1

        N_gamma = gamma.shape[0]
        assignments = []

        for rc in read_classifications:
            if rc.read_id not in gamma_map:
                assignments.append(ReadAssignment(read_id=rc.read_id))
                continue

            idx = gamma_map[rc.read_id]
            if idx >= N_gamma:
                assignments.append(ReadAssignment(read_id=rc.read_id))
                continue

            g_row = gamma[idx]
            rd = read_data[idx]

            posteriors = {}
            for k in range(K):
                if g_row[k] > 1e-6:
                    posteriors[strain_ids[k]] = float(g_row[k])

            # Determine gene family from best candidate
            best_gene = None
            if rd.strain_indices:
                candidate_probs = [g_row[k] for k in rd.strain_indices]
                best_cand_idx = int(np.argmax(candidate_probs))
                if best_cand_idx < len(rd.gene_indices):
                    g = rd.gene_indices[best_cand_idx]
                    if 0 <= g < len(gene_family_ids):
                        best_gene = gene_family_ids[g]

            is_mobile = bool(
                K_total > K and K > 0 and g_row[K] > np.max(g_row[:K])
            )

            assignments.append(
                ReadAssignment(
                    read_id=rc.read_id,
                    posteriors=posteriors,
                    gene_family_id=best_gene,
                    is_mobile=is_mobile,
                )
            )

        return assignments

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _build_phi_priors(
        self, K: int, G: int, reference_phi: np.ndarray | None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build Beta prior hyperparameter matrices.

        If reference annotations are available, use them to set stronger
        priors (a=2.0, b=1.01 for present genes; a=1.01, b=2.0 for absent).
        Otherwise use weak uniform priors.
        """
        cfg = self.config
        prior_a = np.full((K, G), cfg.beta_a)
        prior_b = np.full((K, G), cfg.beta_b)

        if reference_phi is not None:
            # Strengthen priors based on reference annotations
            present = reference_phi > 0.5
            prior_a[present] = max(cfg.beta_a, 2.0)
            prior_b[present] = 1.01
            absent = reference_phi <= 0.5
            prior_a[absent] = 1.01
            prior_b[absent] = max(cfg.beta_b, 2.0)

        return prior_a, prior_b

    def _empty_profile(
        self,
        strain_ids: list[str],
        gene_family_ids: list[str],
        n_reads_total: int,
    ) -> StrainProfile:
        """Return an empty StrainProfile when no reads are classified."""
        K = len(strain_ids)
        G = len(gene_family_ids)
        K_total = K + 1 if self.config.mobile_pool else K

        pi = np.ones(K_total) / K_total
        phi = np.full((K, G), 0.5)

        return StrainProfile(
            species_id="unknown",
            strain_ids=list(strain_ids),
            pi=pi,
            pi_ci_lower=np.zeros(K_total),
            pi_ci_upper=np.ones(K_total),
            phi=phi,
            phi_ci_lower=np.zeros((K, G)),
            phi_ci_upper=np.ones((K, G)),
            gene_family_ids=list(gene_family_ids),
            mobile_pool_proportion=1.0 / K_total if self.config.mobile_pool else 0.0,
            fisher_information=np.zeros((K, K)),
            convergence=ConvergenceDiagnostics(),
            bic=np.inf,
            n_reads_classified=0,
            n_reads_total=n_reads_total,
        )
