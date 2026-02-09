#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "em.h"
#include "utils.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "  FAIL: %s (line %d)\n", msg, __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "  FAIL: %s: %.6f != %.6f (tol=%.6f) (line %d)\n", \
                msg, (a), (b), (tol), __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

/* Helper: generate synthetic EM reads for a 2-species mixture */
static em_read_t *make_reads_2species(int n_reads, int n_markers,
                                       double w0, double w1,
                                       double *bias, /* [2*M] */
                                       uint64_t seed, int *out_n) {
    hs_rng_t rng;
    hs_rng_seed(&rng, seed);

    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_reads, sizeof(em_read_t));
    int count = 0;

    for (int r = 0; r < n_reads; r++) {
        int m = r % n_markers;
        /* Compute probabilities */
        double p0 = w0 * (bias ? bias[0 * n_markers + m] : 1.0);
        double p1 = w1 * (bias ? bias[1 * n_markers + m] : 1.0);
        double total = p0 + p1;
        p0 /= total; p1 /= total;

        /* Assign read */
        double u = hs_rng_uniform(&rng);
        int true_sp = (u < p0) ? 0 : 1;

        reads[count].marker_idx = m;
        reads[count].n_candidates = 2;
        reads[count].species_indices = (int *)hs_malloc(2 * sizeof(int));
        reads[count].containments = (double *)hs_malloc(2 * sizeof(double));
        reads[count].species_indices[0] = 0;
        reads[count].species_indices[1] = 1;
        /* Simulate containment: high for true species, low for false */
        if (true_sp == 0) {
            reads[count].containments[0] = 0.8 + hs_rng_uniform(&rng) * 0.2;
            reads[count].containments[1] = hs_rng_uniform(&rng) * 0.2;
        } else {
            reads[count].containments[0] = hs_rng_uniform(&rng) * 0.2;
            reads[count].containments[1] = 0.8 + hs_rng_uniform(&rng) * 0.2;
        }
        count++;
    }
    *out_n = count;
    return reads;
}

static void test_em_basic_50_50(void) {
    printf("  test_em_basic_50_50...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 42, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 3;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    ASSERT(res->converged, "Marked as converged");

    /* Both weights should be near 0.5 */
    ASSERT_NEAR(res->w[0], 0.5, 0.1, "w[0] near 0.5");
    ASSERT_NEAR(res->w[1], 0.5, 0.1, "w[1] near 0.5");
    ASSERT(res->w[0] + res->w[1] > 0.99, "Weights sum to ~1");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_90_10(void) {
    printf("  test_em_90_10...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.9, 0.1, NULL, 123, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 3;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* Dominant species should have w > 0.7 */
    double dominant = res->w[0] > res->w[1] ? res->w[0] : res->w[1];
    double minor = res->w[0] < res->w[1] ? res->w[0] : res->w[1];
    ASSERT(dominant > 0.7, "Dominant species w > 0.7");
    ASSERT(minor < 0.3, "Minor species w < 0.3");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_bias_recovery(void) {
    printf("  test_em_bias_recovery...\n");
    /* 50/50 mixture with 3x PCR bias on species 1 at marker 0 */
    double bias[6] = { 1.0, 1.0, 1.0,   /* species 0: no bias */
                       3.0, 1.0, 1.0 };  /* species 1: 3x at COI */

    int n_reads;
    em_read_t *reads = make_reads_2species(6000, 3, 0.5, 0.5, bias, 456, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 5;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged with bias");

    /* Despite biased reads, EM should recover near 50/50 weights */
    ASSERT_NEAR(res->w[0], 0.5, 0.15, "w[0] near 0.5 despite bias");
    ASSERT_NEAR(res->w[1], 0.5, 0.15, "w[1] near 0.5 despite bias");

    /* Bias should be recovered in b matrix */
    /* Species 1, marker 0 should have higher bias than marker 1 */
    ASSERT(res->b[1 * 3 + 0] > res->b[1 * 3 + 1] * 0.5,
           "Bias pattern partially recovered");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_confidence_intervals(void) {
    printf("  test_em_confidence_intervals...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 789, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* CIs should contain the true value (0.5) */
    ASSERT(res->w_ci_lo[0] < 0.5 && res->w_ci_hi[0] > 0.5,
           "CI[0] contains true value");
    ASSERT(res->w_ci_lo[1] < 0.5 && res->w_ci_hi[1] > 0.5,
           "CI[1] contains true value");

    /* CIs should be reasonable width */
    double width0 = res->w_ci_hi[0] - res->w_ci_lo[0];
    ASSERT(width0 < 0.5, "CI width < 0.5");
    ASSERT(width0 > 0.001, "CI width > 0.001");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_bic(void) {
    printf("  test_em_bic...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(1000, 3, 0.5, 0.5, NULL, 321, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    ASSERT(isfinite(res->bic), "BIC is finite");
    ASSERT(res->log_likelihood < 0, "Log-likelihood is negative");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_single_species(void) {
    printf("  test_em_single_species...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(1000, 3, 1.0, 0.0, NULL, 654, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    /* Dominant species should get almost all weight */
    double dominant = res->w[0] > res->w[1] ? res->w[0] : res->w[1];
    ASSERT(dominant > 0.8, "Single species gets most weight");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_d_factors(void) {
    printf("  test_em_d_factors...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 111, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* DNA yield factors should be near 1.0 for equal species */
    ASSERT(res->d[0] > 0.3 && res->d[0] < 3.0, "d[0] reasonable range");
    ASSERT(res->d[1] > 0.3 && res->d[1] < 3.0, "d[1] reasonable range");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_single_marker_mode(void) {
    printf("  test_em_single_marker_mode...\n");
    /* Simulate a single-marker (16S) scenario:
     * Species 0 has 2x mito CN vs species 1, so DNA yield is biased.
     * True composition is 50/50, but without mito CN correction the EM
     * would over-estimate the high-CN species from raw read counts. */
    int n_total = 3000;
    int n_reads = 0;
    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_total, sizeof(em_read_t));

    hs_rng_t rng;
    hs_rng_seed(&rng, 999);

    /* Species 0: mito_cn=2000, species 1: mito_cn=1000
     * With 50/50 true weight, species 0 produces 2x more reads */
    double mito_cn[2] = { 2000.0, 1000.0 };
    double eff0 = 0.5 * mito_cn[0];
    double eff1 = 0.5 * mito_cn[1];
    double total_eff = eff0 + eff1;
    double p0 = eff0 / total_eff;  /* ~0.667 of reads from sp0 */

    for (int r = 0; r < n_total; r++) {
        double u = hs_rng_uniform(&rng);
        int true_sp = (u < p0) ? 0 : 1;

        reads[n_reads].marker_idx = 0;  /* single marker */
        reads[n_reads].n_candidates = 2;
        reads[n_reads].species_indices = (int *)hs_malloc(2 * sizeof(int));
        reads[n_reads].containments = (double *)hs_malloc(2 * sizeof(double));
        reads[n_reads].species_indices[0] = 0;
        reads[n_reads].species_indices[1] = 1;
        if (true_sp == 0) {
            reads[n_reads].containments[0] = 0.85 + hs_rng_uniform(&rng) * 0.15;
            reads[n_reads].containments[1] = hs_rng_uniform(&rng) * 0.15;
        } else {
            reads[n_reads].containments[0] = hs_rng_uniform(&rng) * 0.15;
            reads[n_reads].containments[1] = 0.85 + hs_rng_uniform(&rng) * 0.15;
        }
        n_reads++;
    }

    /* Without mito CN correction: d is free, should converge near raw proportions */
    em_config_t cfg_no_cn = em_config_default();
    cfg_no_cn.mito_copy_numbers = NULL;
    int amp_lens[2] = { 460, 460 };

    em_result_t *res_no = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg_no_cn);
    ASSERT(res_no != NULL, "EM (no CN) converged");
    /* Without correction, species 0 gets inflated weight (~0.6+) */
    double w0_no = res_no->w[0] > res_no->w[1] ? res_no->w[0] : res_no->w[1];
    ASSERT(w0_no > 0.55, "Without mito CN: dominant species overestimated");

    /* With mito CN correction: d fixed from known CN values */
    em_config_t cfg_cn = em_config_default();
    cfg_cn.mito_copy_numbers = mito_cn;

    em_result_t *res_cn = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg_cn);
    ASSERT(res_cn != NULL, "EM (with CN) converged");
    /* With correction, both weights should be closer to 0.5 */
    ASSERT_NEAR(res_cn->w[0], 0.5, 0.15, "With mito CN: w[0] near 0.5");
    ASSERT_NEAR(res_cn->w[1], 0.5, 0.15, "With mito CN: w[1] near 0.5");

    /* d should reflect mito CN ratio (2000/1000 -> d[0]/d[1] ~ 2) */
    double d_ratio = res_cn->d[0] / res_cn->d[1];
    ASSERT(d_ratio > 1.5 && d_ratio < 2.5,
           "d ratio reflects mito CN ratio (~2)");

    /* b should be 1.0 in single-marker mode */
    ASSERT_NEAR(res_cn->b[0], 1.0, 0.01, "b[0]=1.0 in single-marker mode");
    ASSERT_NEAR(res_cn->b[1], 1.0, 0.01, "b[1]=1.0 in single-marker mode");

    em_result_destroy(res_no);
    em_result_destroy(res_cn);
    em_reads_free(reads, n_reads);
}

int main(void) {
    printf("=== test_em ===\n");
    test_em_basic_50_50();
    test_em_90_10();
    test_em_bias_recovery();
    test_em_confidence_intervals();
    test_em_bic();
    test_em_single_species();
    test_em_d_factors();
    test_em_single_marker_mode();
    printf("=== %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
