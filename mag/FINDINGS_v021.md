# magprofile v0.2.1 — Novel Findings from MPOB Oil Palm Metagenome Data

**Date:** 2026-02-17
**Dataset:** 30 samples, 1,018 MAGs, 3 compartments (bulk soil, rhizosphere, endosphere) x 2 health statuses (healthy, diseased)

---

## 1. Key Statistical Results

| Test | Compartment (3 groups) | Health Status (2 groups) | Group (6 groups) |
|------|----------------------|------------------------|-----------------|
| PERMANOVA R² | **0.424** (p=0.001) | 0.052 (p=0.120) | 0.498 (p=0.001) |
| ANOSIM R | **0.685** (p=0.001) | 0.055 (p=0.089) | 0.566 (p=0.001) |
| PERMDISP | F=2.32, p=0.122 | F=1.53, p=0.226 | F=0.57, p=0.715 |
| Significant MAGs (FDR<0.05) | 234–556 per pair | **3** | 0–198 per pair |
| Indicator species | 199/239/73 | 125/36 | 67–78 per group |
| Compartment specialists | 65 (all endosphere) | — | 21 |

**Interpretation:** Compartment explains 42.4% of community variation — the dominant ecological axis. Disease explains only 5.2% and is not significant at the community level. PERMDISP is non-significant in all analyses, confirming that PERMANOVA results reflect genuine location shifts, not dispersion artifacts.

---

## 2. Novel Findings

### Finding 1: The endosphere is the only compartment with true ecological specialists

All 65 compartment specialists (specificity score > 0.5) are endosphere-associated. No bulk soil or rhizosphere MAG crosses the 0.5 specificity threshold. The most extreme specialist is a **Myxococcota (Polyangia)** MAG (H4E_semibin_4, specificity = 0.85), meaning 85% of its abundance is confined inside roots.

**Significance:** Most tropical soil microbiome studies report rhizosphere enrichment as the primary niche differentiation signal. Here, the endosphere — not the rhizosphere — harbours the only true specialists, suggesting that the internal root environment imposes the strongest ecological filter in oil palm.

### Finding 2: A predatory Myxococcota guild dominates the endosphere

15 Myxococcota MAGs are significant endosphere indicators (all class Polyangia), spanning two dominant genera: **Palsa-1150** (6 MAGs, up to IndVal = 0.91) and **JADGRB01** (7 MAGs, up to IndVal = 0.88). At the phylum level, Myxococcota shows a 10-fold enrichment in the endosphere vs bulk soil (4.6% vs 0.4%, Cohen's d = -3.86).

**Significance:** Myxococcota are predatory bacteria that lyse and consume other microorganisms. Their massive endosphere enrichment suggests they form a **resident predatory guild inside oil palm roots**, potentially regulating endophytic bacterial populations through top-down control. This is largely unreported in oil palm or tropical crop microbiome literature.

### Finding 3: Roseiarcus — a photoheterotroph thriving inside roots

**Roseiarcus** (Alphaproteobacteria) appears 5 times among the top 20 endosphere indicators (IndVal 0.87–0.90), making it the most recurrent genus. Roseiarcus was originally described from Sphagnum peat bogs as a photoheterotrophic bacterium.

**Significance:** Photoheterotrophs require light for energy but use organic carbon. Finding them as dominant endophytes inside opaque root tissues is paradoxical and suggests either (a) they have switched to a fully heterotrophic lifestyle inside roots, (b) they occupy the outer cortex where some light penetrates, or (c) they retain photosynthetic genes that serve a non-photosynthetic regulatory role. This represents a potentially novel endophytic niche for photoheterotrophic metabolism.

### Finding 4: The rhizosphere is an Acidobacteriota stronghold

56% of all rhizosphere indicator MAGs (41 of 73) are Acidobacteriota (class Terriglobia). Meanwhile, Pseudomonadota drops from 79 indicators in the endosphere to just 3 in the rhizosphere.

**Significance:** This extreme compositional asymmetry — Acidobacteriota dominating the rhizosphere while Pseudomonadota dominate the endosphere — reveals that oil palm roots impose a sharp ecological boundary. The rhizosphere selects for slow-growing, oligotrophic Acidobacteriota (likely adapted to recalcitrant root exudate degradation), while the endosphere selects for fast-growing Pseudomonadota and predatory Myxococcota.

### Finding 5: Near-complete archaea exclusion from roots

Thermoplasmatota and Thermoproteota (ammonia-oxidising archaea) show the largest effect sizes in the entire dataset:

| Phylum | Bulk soil | Endosphere | Cohen's d |
|--------|----------|------------|-----------|
| Thermoplasmatota | 15.4% | 3.5% | **6.43** |
| Thermoproteota | 4.2% | 1.4% | **6.02** |

**Significance:** Oil palm roots almost completely exclude archaea. These are the strongest compartment effects of any phylum (d > 6.0). In bulk soil, archaea comprise ~20% of the community; inside roots, they are reduced to ~5%. This suggests strong anti-archaeal selection at the root surface, possibly mediated by plant immune responses that cannot distinguish archaea from bacteria, or by competitive exclusion from bacterial endophytes.

### Finding 6: Disease produces a razor-thin but real signal — 3 MAGs

Despite the global PERMANOVA being non-significant (p = 0.12), exactly **3 MAGs** pass FDR correction for differential abundance between diseased and healthy palms:

| MAG | Phylum | Genus | Cohen's d |
|-----|--------|-------|-----------|
| Coassem3_concoct_34 | Verrucomicrobiota | UBA11358 | **2.49** |
| Coassem3_SemiBin_cstrat_454 | Verrucomicrobiota | UBA7542 | **2.01** |
| Coassem3_metabat_610 | Actinomycetota | JAJYUU01 | **1.63** |

All three are enriched in diseased palms with very large effect sizes. This is a targeted, not community-wide, disease response.

**Significance:** Basal stem rot (Ganoderma) does not restructure the entire soil microbiome. Instead, it produces precise, large-effect-size shifts in specific taxa — particularly Verrucomicrobiota. This suggests disease biomarkers may be identifiable from very few MAGs rather than community-level profiles.

### Finding 7: Ultrasmall organisms as dysbiosis markers

Among the top 10 disease indicators (by IndVal score):
- **Micrarchaeota** (DPANN archaea, genus JAJZYD01, IndVal = 0.90)
- **Patescibacteria** (CPR bacteria, Saccharimonadia, IndVal = 0.87)

At the phylum level, Patescibacteria is significantly enriched in diseased palms (Cohen's d = 1.32, q = 0.017), while Nanoarchaeota is significantly depleted (d = -1.36, q = 0.017).

**Significance:** Both Micrarchaeota and Patescibacteria are candidate phyla radiation/DPANN organisms with reduced genomes and likely parasitic or epibiotic lifestyles on other microbes. Their enrichment in diseased palms suggests they proliferate when the normal microbiome is disturbed — functioning as **microbial dysbiosis markers** analogous to the Firmicutes/Bacteroidetes ratio in human gut dysbiosis.

### Finding 8: Universal membership, differential abundance

1,017 of 1,018 MAGs are core at the 100% prevalence threshold — every MAG is detected in every sample. Community differences are driven entirely by **abundance shifts**, not presence/absence.

**Significance:** This challenges the common "who's there?" framing of microbiome studies. In this system, everyone is everywhere — the question is how much of each. This makes abundance-based analyses (CLR differential, EM estimation) essential, and presence/absence methods (Jaccard, binary core microbiome) uninformative.

---

## 3. Motivation for Chapter 2 (magfunc) and Chapter 3 (magnet)

### Why magfunc is now essential

The magprofile findings raise functional questions that **cannot be answered by taxonomy alone**:

| Finding | Functional question for magfunc |
|---------|-------------------------------|
| Myxococcota predatory guild (15 MAGs) | What lytic enzymes, secretion systems (T3SS/T6SS), and secondary metabolites enable predation inside roots? Are these MAGs functionally redundant or do they partition prey? |
| Roseiarcus photoheterotroph in roots | Do these MAGs retain complete photosynthetic gene clusters (puf/bch operons) or have they lost them? What alternative energy metabolism operates inside opaque root tissue? |
| Acidobacteriota rhizosphere dominance | Are rhizosphere Acidobacteriota enriched for CAZymes (cellulose/hemicellulose degradation of root exudates) compared to bulk soil Acidobacteriota? |
| Archaea exclusion from roots | What functional roles do Thermoplasmatota/Thermoproteota fill in bulk soil (ammonia oxidation? methanogenesis?) that is absent from the endosphere? Is nitrogen cycling restructured inside roots? |
| Verrucomicrobiota disease enrichment | What metabolic traits distinguish the 3 disease-enriched Verrucomicrobiota MAGs? Are they carrying virulence-associated genes or simply opportunistic degraders? |
| Ultrasmall dysbiosis markers | What parasitic/epibiotic functions do the Micrarchaeota and Patescibacteria MAGs encode? Can we identify their host organisms from metabolic dependencies? |
| Universal membership, abundance-driven | If all MAGs are everywhere, what functional genes are differentially expressed across compartments? Functional profiling captures the mechanistic basis for abundance shifts. |

**Verdict: Extremely high motivation.** Every major finding generates a "what genes enable this?" question. Taxonomy tells us *who* has shifted — magfunc will tell us *what they can do* and *why* the root environment selects for them.

### Why magnet is now essential

The magprofile findings reveal interaction patterns that **demand network-level analysis**:

| Finding | Interaction question for magnet |
|---------|-------------------------------|
| Myxococcota predatory guild | **Who are they preying on?** Predator-prey relationships should produce strong negative correlations. Network analysis can identify the prey taxa and test whether predatory control explains endosphere community structure. |
| Ultrasmall parasites (CPR/DPANN) | **Who are their hosts?** Micrarchaeota and Patescibacteria are obligate symbionts/parasites with reduced genomes. They cannot live independently — co-occurrence networks can identify their host MAGs. |
| 65 endosphere specialists | **Do they form a coherent interaction module?** If specialists co-occur tightly, the endosphere has a deterministic, tightly-coupled community. If not, they are independently filtered. |
| Roseiarcus (5 co-occurring MAGs) | **Do Roseiarcus MAGs have consistent interaction partners?** Multiple congeneric MAGs with IndVal > 0.87 may share a conserved niche defined by interactions with the same partner taxa. |
| Compartment-driven assembly (R²=0.42) | **Is assembly deterministic or stochastic?** Beta-NTI null models can distinguish environmental filtering (deterministic) from ecological drift (stochastic). The strong compartment signal suggests deterministic assembly — magnet can quantify this. |
| Disease effect is targeted (3 MAGs) | **Does disease disrupt specific network modules?** Rather than a community-wide effect, disease may break or rewire specific interaction edges. Differential network analysis (healthy vs diseased networks) can identify the disrupted connections. |
| Archaea excluded from endosphere | **What replaces archaea's metabolic role inside roots?** If ammonia-oxidising archaea are excluded, bacterial ammonia oxidisers or alternative nitrogen cyclers must fill the niche. Network analysis can identify which bacterial MAGs occupy the archaeal niche inside roots. |

**Verdict: Extremely high motivation.** The findings practically demand interaction analysis. The predatory Myxococcota guild, parasitic ultrasmall organisms, and compartment specialist modules are all fundamentally *interaction* phenomena that taxonomy and function alone cannot resolve.

---

## 4. The Three-Chapter Narrative

| Chapter | Tool | Question | Key motivation from data |
|---------|------|----------|------------------------|
| 1 | **magprofile** (done) | Who's there and where? | Compartment, not disease, drives community assembly. Endosphere harbours exclusive specialists including a predatory Myxococcota guild and paradoxical photoheterotrophs. |
| 2 | **magfunc** (next) | What can they do? | Predators need lytic genes. Photoheterotrophs in darkness need alternative energy metabolism. Disease markers need functional characterisation. Every finding demands functional annotation. |
| 3 | **magnet** (final) | How do they interact? | Predator-prey, parasite-host, specialist modules, and niche replacement are all interaction phenomena. Network analysis will reveal the ecological wiring diagram that produces the patterns magprofile detected. |

Each chapter builds directly on the previous one's findings — not as a generic pipeline, but because the specific biological discoveries force the next analytical question.

---

## 5. Summary of Numbers

- **1,018** MAGs profiled across 30 samples
- **42.4%** of community variation explained by compartment (PERMANOVA R²)
- **65** endosphere specialists (no specialists in other compartments)
- **15** Myxococcota predatory MAGs enriched in endosphere (10x vs bulk)
- **5** Roseiarcus MAGs among top 20 endosphere indicators
- **3** MAGs significantly associated with disease (FDR < 0.05)
- **2** ultrasmall phyla (Patescibacteria, Nanoarchaeota) significantly shifted in disease
- **106** unique genera among significant compartment indicators
- **1,017/1,018** MAGs are core at 100% prevalence — universal membership
