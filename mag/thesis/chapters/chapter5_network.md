# Chapter 5: Co-occurrence Networks and Keystone Taxa in Oil Palm Soils

**Estimated length:** 25–35 pages | **Writing priority:** 5th (June 2026)

**Relationship to Paper 2:** This chapter expands the network analysis component of Paper 2 (Pertanika JTAS) into thesis format. Additional content includes: expanded methodology, detailed module-by-module analysis, extended keystone characterisation, and thesis-specific framing connecting to Objective 3.

**Source:** Paper 2 manuscript (`thesis/papers/paper2_functional_network.md`) + `v043_analysis_report.txt` Sections 3, 6, 7

---

## 5.1 Introduction

[Adapted from Paper 2 Introduction, focused on network analysis:]
- Explicit link to Objective 3
- Why network analysis is necessary given Chapter 3 and 4 findings
- Motivation table: each ecological/functional finding → network question (from FINDINGS_v021.md Section 3)
- Overview of phi proportionality, modularity, keystone, hub-bridge approaches

---

## 5.2 Materials and Methods

### 5.2.1 Network Construction

[Same as Paper 2, expanded with:]
- CLR transformation details
- Phi proportionality computation: formula, implementation
- Prevalence filtering (50% minimum)
- Threshold selection rationale

### 5.2.2 Threshold Sensitivity Analysis

- Five percentiles tested: 1st, 3rd, 5th, 10th, 20th
- Metrics evaluated: edges, density, modularity, modules
- Selection criteria: balance resolution with edge count

### 5.2.3 Module Detection and Null Model Validation

[Same as Paper 2, expanded with:]
- Louvain algorithm: resolution parameter
- Configuration model null: degree sequence preservation
- 1,000 random networks
- Z-score computation

### 5.2.4 Keystone Identification

[Same as Paper 2, expanded with:]
- Composite scoring formula
- Centrality metric definitions (betweenness, closeness, degree)
- Abundance inverse-weighting
- Threshold for keystone classification

### 5.2.5 Hub-Bridge Classification

[Same as Paper 2, expanded with:]
- z-P computation per MAG
- Classification thresholds and justification
- Module membership assignment

### 5.2.6 Compartment-Specific and Health-Specific Networks

[Same as Paper 2, expanded with:]
- Sample subsetting strategy
- Recalculation of phi values per subset
- Differential edge identification algorithm

### 5.2.7 Software

[magprofile v0.43, magnet module]

---

## 5.3 Results

### 5.3.1 Global Network Properties

**Table 5.1.** Global network summary:

| Property | Value |
|----------|------:|
| MAGs | 1,018 |
| Edges | 25,883 |
| Density | 0.050 |
| Modularity (*Q*) | 0.226 |
| Modules | 258 |
| Phi threshold | 0.0532 (5th percentile) |

### 5.3.2 Threshold Sensitivity

**Table 5.2.** Network properties across thresholds:

| Percentile | Threshold | Edges | Modularity | Modules |
|-----------:|----------:|------:|-----------:|--------:|
| 1% | 0.029 | 5,177 | 0.306 | 521 |
| 3% | 0.043 | 15,530 | 0.250 | 341 |
| 5% | 0.053 | 25,883 | 0.226 | 258 |
| 10% | 0.073 | 51,766 | 0.189 | 163 |
| 20% | 0.106 | 103,531 | 0.158 | 78 |

Figure 5.1: Modularity vs threshold curve showing robustness

### 5.3.3 Null Model Validation

- Observed Q = 0.226
- Null mean = 0.073 ± 0.002
- **z = 81.36, p < 0.001**
- Figure 5.2: Null distribution with observed value indicated

Context: z = 81.36 exceeds values typically reported (z = 5–30; Li et al., 2022). The extreme value reflects large sample size (1,018 MAGs) giving precise null estimates (SD = 0.002).

### 5.3.4 Module Composition

**Table 5.3.** Composition of top modules:

| Module | MAGs | Dominant phylum | Fraction | Key features |
|-------:|-----:|-----------------|:--------:|:-------------|
| 0 | 341 | Acidobacteriota | 33.4% | Keystone hub, 8 Myxococcota |
| 1 | 275 | Pseudomonadota | 24.7% | Copiotrophic cluster |
| 2 | 117 | Pseudomonadota | 23.9% | Secondary copiotrophic |
| 3 | 9 | Actinomycetota | 55.6% | Small specialist module |
| 4 | 6 | Pseudomonadota | 66.7% | Small specialist module |

- 258 total modules; 245 are singletons or small (1–2 MAGs)
- 36 phyla represented across all modules
- Module taxonomic coherence analysis: Figure 5.3

### 5.3.5 Keystone Taxa

**Table 5.4.** Top keystone taxa:

| MAG | Phylum | Class | Genus | Score | Betweenness | Abundance (%) |
|-----|--------|-------|-------|------:|------------:|--------------:|
| Coassem3_SemiBin_cstrat_619 | Acidobacteriota | Terriglobia | *Terracidiphilus* | 0.981 | 0.012 | 0.30 |
| Coassem2_maxbin2_438 | Acidobacteriota | Terriglobia | — | 0.896 | — | — |
| Coassem2_maxbin2_192 | Acidobacteriota | Terriglobia | — | 0.883 | — | — |
| Coassem2_maxbin2_300 | Verrucomicrobiota | — | — | 0.848 | — | — |
| Coassem3_SemiBin_cstrat_833 | Chlamydiota | — | JABDCS01 | 0.723 | — | — |

- *Terracidiphilus*: "rare but important" pattern (low abundance, highest centrality)
- Acidobacteriota dominate top keystone positions (3 of top 5)
- Figure 5.4: Keystone score distribution with threshold

### 5.3.6 Myxococcota in the Network

**Table 5.5.** All 26 Myxococcota MAGs in network context:

- 26 total (2.6% of 1,018)
- Dominant class: Polyangia (24/26 = 92%)
- Most common genus: *Palsa-1150* (13 MAGs)
- Keystone Myxococcota: 5 (19%)
- Bridge Myxococcota: 3 (12%)
- In Module 0: 8 MAGs

**Table 5.6.** Bridge Myxococcota:

| MAG | Genus | P | z | Keystone? |
|-----|-------|--:|--:|:---------:|
| Coassem3_metabat_193 | Polyangia | 0.629 | — | Yes |
| Coassem3_SemiBin_cstrat_536 | *Palsa-1150* | 0.643 | — | Yes |
| Coassem3_SemiBin_cstrat_790 | *Palsa-1150* | 0.625 | — | Yes |

- All three bridges are also keystones
- Figure 5.5: Subnetwork of Module 0 highlighting Acidobacteriota-Myxococcota co-occurrence

### 5.3.7 Hub-Bridge Classification

[Expanded from Paper 2]
- z-P scatter plot: Figure 5.6
- Hub nodes: high within-module degree (z > 2.5)
- Bridge nodes: high participation coefficient (P > 0.62)
- Peripheral: majority of MAGs
- Classification counts per phylum

### 5.3.8 Compartment-Specific Networks

**Table 5.7.** Compartment network properties:

| Network | Edges | Density |
|---------|------:|--------:|
| Global | 25,883 | 0.050 |
| Bulk | 141,465 | 0.273 |
| Endosphere | 119,327 | 0.231 |
| Rhizosphere | 117,477 | 0.227 |

- Density gradient: bulk > endosphere > rhizosphere
- Figure 5.7: Side-by-side network visualisations

### 5.3.9 Health-Specific Networks and Structural Dysbiosis

**Table 5.8.** Health-specific network properties:

| Network | Edges | Density |
|---------|------:|--------:|
| Diseased | 60,351 | 0.117 |
| Healthy | 47,416 | 0.092 |
| **Difference** | **+27%** | |

**Table 5.9.** Group-level network properties:

| Group | Edges | Density |
|-------|------:|--------:|
| Diseased bulk | 186,809 | 0.361 |
| Diseased endosphere | 167,729 | 0.324 |
| Diseased rhizosphere | 178,923 | 0.346 |
| Healthy bulk | 187,003 | 0.361 |
| Healthy endosphere | 151,894 | 0.293 |
| Healthy rhizosphere | 159,151 | 0.307 |

- Diseased consistently denser across all compartments (except bulk ≈ equal)
- Edge gain analysis: 90,544 edges in differential network

### 5.3.10 Differential Network Analysis

[Expanded from Paper 2]
- Compartment differential edges:
  - Bulk vs Endosphere: 232,388 edges analysed
  - Bulk vs Rhizosphere: 193,986 edges
  - Endosphere vs Rhizosphere: 186,843 edges
- Health differential edges:
  - Diseased vs Healthy: 90,544 edges (conserved + gained + lost)
- Figure 5.8: Differential edge visualisation
- High conservation between endosphere and rhizosphere → shared root-associated networks

---

## 5.4 Discussion

### 5.4.1 Non-Random Modular Architecture

[Expanded from Paper 2]
- z = 81.36 as strongest validation in soil microbiome literature
- Comparison with other z-scores: Li et al. (2022) review
- Ecological interpretation: deterministic assembly, not ecological drift
- Module size distribution: power law or truncated?

### 5.4.2 The Acidobacteriota-Myxococcota Trophic Cascade

[Major expansion from Paper 2]
- Acidobacteriota as keystone: Module 0, 341 MAGs, *Terracidiphilus* (0.981)
- *Terracidiphilus* characterisation: cellulose/chitin degradation (García-Fraile et al., 2016)
- Myxococcota as bridges: 3 Polyangia MAGs connecting modules
- Predator-prey dynamics: Zhou et al. (2020) negative correlation
- Nutrient recycling model:
  - Acidobacteriota degrade complex C → release simple sugars
  - Myxococcota prey on Acidobacteriota → release N, P, C back to soil
  - This "microbial loop" drives nutrient cycling
- Comparison with marine viral shunt
- Dai et al. (2023): Myxococcota diversity linked to multi-nutrient cycling

### 5.4.3 Structural Dysbiosis under Disease

[Major expansion from Paper 2]
- 27% denser network in diseased soils
- No compositional change (Chapter 3) + no functional change (Chapter 4) + structural change = structural dysbiosis
- Hypotheses:
  - (a) Loss of niche partitioning → homogenisation of co-occurrence
  - (b) Cooperative stress response → more cross-talk
- Literature comparison:
  - Gao et al. (2021): Fusarium increased complexity (supports)
  - Shi et al. (2021): Verticillium increased complexity (supports)
  - Deng et al. (2021): Banana Fusarium decreased complexity (contrasts)
- Implications: network density as early-warning indicator for BSR

### 5.4.4 Myxococcota and Disease Response

[Expanded from Paper 2]
- Expanded interaction profiles in diseased networks
- Bacteroidota-Myxococcota and Verrucomicrobiota-Myxococcota links
- Biocontrol literature: Ye et al. (2020), Zhao et al. (2024), Kuang et al. (2023)
- Predatory recruitment hypothesis

### 5.4.5 Network Analysis Complements Community Profiling

- Standard approaches (PERMANOVA, differential abundance) miss structural dysbiosis
- Network analysis detects disease signal invisible to compositional methods
- Recommendation: network analysis as standard component of disease microbiome studies

### 5.4.6 Limitations

[Same as Paper 2, expanded:]
- Co-occurrence ≠ interaction
- Correlation does not imply causation
- Threshold sensitivity
- Small sample size for group-level networks (n = 5 per group)
- Need for experimental validation (predation assays, stable isotope probing)

---

## 5.5 Conclusions

[Mapped explicitly to Objective 3:]
1. The oil palm microbiome forms 258 non-random modules (z = 81.36), indicating deterministic assembly
2. *Terracidiphilus* (Acidobacteriota) is the top keystone taxon; Myxococcota serve as inter-module bridges
3. The Acidobacteriota-Myxococcota trophic cascade anchors community metabolism and connectivity
4. Disease produces structural dysbiosis (27% denser) without compositional or functional change
5. These results support the network layer of the insurance hypothesis

---

## 5.6 References

[Superset of Paper 2 network references + additional thesis-specific citations]
