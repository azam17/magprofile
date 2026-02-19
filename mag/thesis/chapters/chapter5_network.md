# Chapter 5: Co-occurrence Networks and Keystone Taxa in Oil Palm Root-Associated Microbiomes

*This chapter addresses Objective 3 and corresponds to the network analysis component of Paper 2.*

---

## 5.1 Introduction

[Expanded from Paper 2 Introduction, focused specifically on network questions. ~3-4 pages]

Key points:
- Chapters 3-4 established compartment-driven composition and function; this chapter examines interaction architecture
- Table linking ecological findings to interaction questions (from FINDINGS_v021.md Section 3):
  - Myxococcota guild → who are they preying on?
  - Ultrasmall parasites → who are their hosts?
  - 65 endosphere specialists → do they form a coherent interaction module?
  - Roseiarcus (5 MAGs) → consistent interaction partners?
  - Compartment dominance → deterministic or stochastic assembly?
  - Disease effect (3 MAGs) → specific network module disruption?
- Co-occurrence network analysis background
- Specific hypotheses:
  - H1: Network structure is significantly non-random (null model test)
  - H2: Modules reflect taxonomic coherence
  - H3: Keystone taxa are low-abundance but structurally critical
  - H4: Disease alters network topology without changing composition

## 5.2 Materials and Methods

### 5.2.1 Network construction

**Source:** Paper 2 Methods

- 1,018 MAGs (both HQ and MQ) used for abundance-based network
- Phi proportionality (Lovell et al., 2015) on CLR-transformed abundances
- Prevalence filter: ≥50% of samples
- Threshold: 5th percentile of phi distribution (|phi| < 0.0532)
- Threshold sensitivity: 1st, 3rd, 5th, 10th, 20th percentiles tested

### 5.2.2 Module detection

- Louvain algorithm (Blondel et al., 2008)
- Modularity Q (Newman, 2006)
- Module composition analysis (phylum proportions)

### 5.2.3 Null model validation

- 1,000 random networks via configuration model (preserved degree sequence)
- z-score = (Q_observed - Q_null_mean) / Q_null_SD
- p-value from null distribution

### 5.2.4 Keystone taxa identification

- Composite score from normalised betweenness, closeness, and degree centrality
- Inverse-weighted by mean relative abundance
- "Rare but structurally important" criterion

### 5.2.5 Hub-bridge classification

- Guimerà-Amaral (2005) z-P framework
- Within-module degree z-score (hub threshold: z > 2.5)
- Participation coefficient P (bridge threshold: P > 0.62)
- Four-quadrant classification: hub, connector hub, bridge, peripheral

### 5.2.6 Compartment-specific and health-specific networks

- Separate networks for bulk, endosphere, rhizosphere (subset samples, recalculate phi)
- Separate networks for diseased, healthy
- Network density comparison

### 5.2.7 Differential network analysis

- Edge set comparison: conserved, gained, lost edges between health states
- Compartment differential edges

## 5.3 Results

### 5.3.1 Global network properties

**Source:** v043_analysis_report.txt Section 3.1

- 1,018 MAGs, 25,883 edges
- Network density: 0.050
- Modularity Q: 0.226
- 258 modules
- Phi threshold: 0.0532 (5th percentile)

### 5.3.2 Non-random modular structure (null model)

**Source:** v043_analysis_report.txt Section 3.2

- Observed modularity: 0.226
- Null model mean: 0.073 ± 0.002
- z-score: 81.36
- p < 0.001

81 standard deviations above random — the community structure is overwhelmingly non-random.

Context: Li et al. (2022) define z > 2 as significant, modularity > 0.4 as modular. z = 81 far exceeds typical reports (z = 5-30), reflecting the statistical power of 1,018 MAGs.

### 5.3.3 Threshold sensitivity

**Source:** v043_analysis_report.txt Section 3.3

| Percentile | Threshold | Edges | Modularity | Modules |
|----------:|----------:|------:|-----------:|--------:|
| 1% | 0.029 | 5,177 | 0.306 | 521 |
| 3% | 0.043 | 15,530 | 0.250 | 341 |
| 5% | 0.053 | 25,883 | 0.226 | 258 |
| 10% | 0.073 | 51,766 | 0.189 | 163 |
| 20% | 0.106 | 103,531 | 0.158 | 78 |

Modularity robust across 20-fold edge count range (0.158–0.306).

### 5.3.4 Compartment-specific networks

**Source:** v043_analysis_report.txt Section 3.4

| Network | Edges | Density |
|---------|------:|--------:|
| Global | 25,883 | 0.050 |
| Bulk | 141,465 | 0.273 |
| Endosphere | 119,327 | 0.231 |
| Rhizosphere | 117,477 | 0.227 |

Density gradient: bulk > endosphere > rhizosphere.

### 5.3.5 Disease increases network density by 27%

**Source:** v043_analysis_report.txt Sections 3.5–3.6

| Network | Edges | Density |
|---------|------:|--------:|
| Diseased | 60,351 | 0.117 |
| Healthy | 47,416 | 0.092 |

27% more edges in diseased networks.

Group-level:
| Group | Edges | Density |
|-------|------:|--------:|
| diseased_bulk | 186,809 | 0.361 |
| diseased_endo | 167,729 | 0.324 |
| diseased_rhizo | 178,923 | 0.346 |
| healthy_bulk | 187,003 | 0.361 |
| healthy_endo | 151,894 | 0.293 |
| healthy_rhizo | 159,151 | 0.307 |

Disease effect strongest in endosphere (0.324 vs 0.293) and rhizosphere (0.346 vs 0.307).

### 5.3.6 Module composition

**Source:** v043_analysis_report.txt Section 3.7

| Module | MAGs | Dominant phylum | Fraction |
|-------:|-----:|-----------------|--------:|
| 0 | 341 | Acidobacteriota | 33.4% |
| 1 | 275 | Pseudomonadota | 24.7% |
| 2 | 117 | Pseudomonadota | 23.9% |
| 3 | 9 | Actinomycetota | 55.6% |
| 4 | 6 | Pseudomonadota | 66.7% |

258 total modules; 245 are singletons or small (1-2 MAGs). 36 phyla represented.

### 5.3.7 Keystone taxa

**Source:** v043_analysis_report.txt Section 3.8

Top keystone: **Terracidiphilus** (Acidobacteriota, Terriglobia)
- MAG: Coassem3_SemiBin_cstrat_619
- Score: 0.981
- Betweenness centrality: 0.012
- Mean abundance: 0.30%
- Classic "rare but important" pattern

Other notable keystones:
- Coassem2_maxbin2_438 (Acidobacteriota, Terriglobia): 0.896
- Coassem2_maxbin2_192 (Acidobacteriota, Terriglobia): 0.883
- Coassem2_maxbin2_300 (Verrucomicrobiota): 0.848
- Coassem3_SemiBin_cstrat_833 (Chlamydiota, JABDCS01): 0.723

Acidobacteriota dominate keystone positions.

### 5.3.8 Hub-bridge classification

**Source:** v043_analysis_report.txt Section 3.9

Guimerà-Amaral z-P framework:
- Hub nodes (z > 2.5): high within-module connectivity
- Bridge nodes (P > 0.62): connect multiple modules
- Most MAGs: peripheral (low z, low P)

### 5.3.9 Myxococcota as network bridges

**Source:** v043_analysis_report.txt Section 6

26 Myxococcota MAGs total:
- 92% class Polyangia
- 13 genus Palsa-1150
- 5 identified as keystones (19%)
- 3 bridges (P = 0.625–0.643)
- 8 in Module 0 (Acidobacteriota hub)

Bridge Myxococcota:
| MAG | Genus | P | Keystone |
|-----|-------|--:|:--------:|
| Coassem3_metabat_193 | Polyangia | 0.629 | Yes |
| Coassem3_SemiBin_cstrat_536 | Palsa-1150 | 0.643 | Yes |
| Coassem3_SemiBin_cstrat_790 | Palsa-1150 | 0.625 | Yes |

### 5.3.10 Differential networks

**Source:** v043_analysis_report.txt Section 3.7 (differential)

Compartment comparisons:
- Bulk vs Endosphere: 232,388 edges analysed
- Bulk vs Rhizosphere: 193,986 edges
- Endo vs Rhizosphere: 186,843 edges

Health comparison:
- Diseased vs Healthy: 90,544 edges (conserved + gained + lost)
- Increased density driven by edge gain, not loss

## 5.4 Discussion

### 5.4.1 Non-random modular organisation

[~2 pages]

- z = 81.36 is the strongest statistical validation of non-random structure
- Comparison with Guimerà & Amaral (2005): z > 2 significant
- Li et al. (2022): most studies z = 5-30
- Extreme value reflects 1,018-MAG statistical power
- Biological interpretation: deterministic assembly processes dominate

### 5.4.2 The Acidobacteriota-Myxococcota trophic axis

[~3 pages. From v043_analysis_report.txt Sections 6.1-6.5]

- Module 0 (341 MAGs): Acidobacteriota anchor
- Terracidiphilus: cellulose/chitin degrader (García-Fraile et al., 2016), keystone score 0.981
- 8 Myxococcota within Module 0: predator-prey co-occurrence
- Zhou et al. (2020): r = -0.495 between phyla (negative correlation = predation)
- Trophic cascade model:
  - Acidobacteriota → complex carbon degradation → community metabolism base
  - Myxococcota → predation → nutrient recycling (C, N, P)
  - Bridges → cross-module connection via diverse prey
- Dai et al. (2023): Myxococcota phylogenetic diversity linked to nutrient cycling
- Goncalves et al. (2024): Acidobacteriota carry PGP genes

### 5.4.3 Structural dysbiosis: disease changes wiring, not components

[~3 pages. From v043_analysis_report.txt Sections 7.1-7.3]

- 27% density increase without compositional or functional change
- Competing hypotheses: (a) loss of niche partitioning, (b) cooperative stress response
- Supporting literature:
  - Gao et al. (2021): Fusarium wilt increased network complexity
  - Shi et al. (2021): infected samples had more complex networks
  - Deng et al. (2021): banana Fusarium — contrasting (fewer connections)
- Myxococcota interaction expansion in disease:
  - New Bacteroidota-Myxococcota, Verrucomicrobiota-Myxococcota links in diseased networks
  - Ye et al. (2020): Corallococcus biocontrol against Fusarium (54-80% reduction)
  - Zhao et al. (2024): Myxococcota enriched in Verticillium wilt
  - Possible predatory defence recruitment

### 5.4.4 Practical implications: network density as early warning

[~1 page]

- Standard PERMANOVA misses the disease signal entirely (p = 0.120)
- Network density captures structural reorganisation invisible to compositional analysis
- Potential early-warning indicator for BSR before visual symptoms appear
- Monitoring approach: periodic metagenomic sampling with network density tracking

## 5.5 Summary

- 25,883 edges organised into 258 non-random modules (z = 81.36)
- Terracidiphilus (Acidobacteriota) is the top keystone (score 0.981)
- Myxococcota serve as inter-module bridges connecting communities via predation
- Disease increases density by 27% = structural dysbiosis
- Network analysis captures disease effects invisible to standard approaches

*This completes the three-pillar analysis: ecology (Chapter 3), function (Chapter 4), and network (Chapter 5). Chapter 6 synthesises these findings into the insurance hypothesis and discusses broader implications.*

---

*[Data sources: v043_analysis_report.txt Sections 3.1–3.9, 6.1–6.5, 7.1–7.3]*
