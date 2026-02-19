# Chapter 4: Functional Capacity and PGPR Potential of Compartment-Stratified MAGs

**Estimated length:** 25–35 pages | **Writing priority:** 4th (May 2026)

**Relationship to Paper 2:** This chapter expands the functional profiling component of Paper 2 (Pertanika JTAS) into thesis format. Additional content includes: expanded methods, deeper results on individual pathways and CAZyme families, extended discussion with trait-by-trait analysis, and thesis-specific framing connecting to Objective 2.

**Source:** Paper 2 manuscript (`thesis/papers/paper2_functional_network.md`) + `v043_analysis_report.txt` Sections 2, 4, 5

---

## 4.1 Introduction

[Adapted from Paper 2 Introduction, focused on functional analysis:]
- Explicit link to Objective 2
- Why functional annotation is necessary given Chapter 3 findings
- Motivation table: each Chapter 3 finding → functional question (from FINDINGS_v021.md Section 3)
- Overview of DRAM, KEGG, CAZyme, PGPR trait approaches

---

## 4.2 Materials and Methods

### 4.2.1 MAG Selection for Functional Annotation

- 567 HQ MAGs (completeness ≥90%, contamination ≤5%) selected for DRAM annotation
- Rationale: gene absence on incomplete genomes is unreliable; MQ MAGs excluded to avoid false negatives
- Tradeoff: reduced statistical power for rare traits (e.g., nif genes: 5–7 carriers out of 567)

### 4.2.2 DRAM Annotation

[Same as Paper 2, with additional detail:]
- DRAM version, database versions (KEGG release, dbCAN version)
- Annotation confidence thresholds
- Quality control: filtering low-confidence annotations

### 4.2.3 Functional Profiling Analyses

#### KO Pathway Abundance
- Per-sample: sum of relative abundances of MAGs carrying each KO
- CLR transformation for differential analysis

#### CAZyme Classification
- dbCAN family assignment
- Six-class breakdown (GH, GT, CE, AA, CBM, PL)
- Per-sample abundance calculation

#### PGPR Trait Screening
- 13-gene panel (Table 4.1)
- Binary presence/absence from KO annotations
- Statistical testing: Kruskal-Wallis with FDR correction

#### Functional Redundancy
- Shannon diversity of function carriers per function per sample
- Comparison across compartments and health states

#### Differential Functional Analysis
- Mann-Whitney U test with FDR correction (q < 0.05)
- All pairwise comparisons: compartment (3), health (1), group (15)

### 4.2.4 Software

[magprofile v0.43, magfunc module]

---

## 4.3 Results

### 4.3.1 Annotation Summary

- 6,233 KO functions detected
- 532 CAZyme families across 6 classes
- 13 PGPR traits screened

**Table 4.2.** CAZyme class distribution (from Paper 2 Table 1):

| Class | Genes | Families | MAGs |
|-------|------:|--------:|-----:|
| GH | 26,069 | 193 | 557 (98.2%) |
| GT | 22,083 | 97 | 567 (100%) |
| CE | 5,155 | 42 | 529 (93.3%) |
| AA | 3,764 | 10 | 530 (93.5%) |
| CBM | 2,061 | 159 | 438 (77.2%) |
| PL | 824 | 31 | 301 (53.1%) |

### 4.3.2 Compartment-Driven Functional Differentiation

[Same as Paper 2, expanded with:]
- Top 20 most differentially abundant KO pathways per comparison: Table 4.3
- Top 20 most differentially abundant CAZyme families: Table 4.4
- Heatmap of pathway abundance across compartments: Figure 4.1
- Heatmap of CAZyme abundance across compartments: Figure 4.2
- Volcano plots for each pairwise comparison: Figure 4.3

Summary counts:
- Bulk vs Endosphere: ~2,000 differential KOs, ~200 CAZymes
- Bulk vs Rhizosphere: ~1,500 KOs, ~150 CAZymes
- Endosphere vs Rhizosphere: ~500 KOs, ~50 CAZymes
- Diseased vs Healthy: 0 KOs, 0 CAZymes

### 4.3.3 PGPR Trait Enrichment

[Same as Paper 2, expanded with:]

**Table 4.5.** Full PGPR trait results (from Paper 2 Table 3)

Detailed analysis per trait:

#### 4.3.3.1 Nitrogen Fixation (*nifHDK*)
- *nifH*: 0 bulk, 5 endo, 1 rhizo (q = 0.031)
- *nifK*: 0, 6, 1 (q = 0.016)
- *nifD*: 0, 4, 1 (q = 0.078, trend)
- Coherent signal across all three nitrogenase subunits
- Identity of nif-carrying MAGs: taxonomic affiliation, compartment profiles

#### 4.3.3.2 ACC Deaminase (*acdS*)
- 5 bulk, 15 endo, 3 rhizo (q = 0.006)
- Strongest PGPR enrichment signal
- Carrier MAG characterisation

#### 4.3.3.3 Phosphate Solubilisation (*gcd*)
- 120 bulk, 117 endo, 84 rhizo (q = 0.006)
- Enriched in bulk soil — free-living P-cycling dominates

#### 4.3.3.4 Other Traits
- *budB*: 162, 131, 103 (q = 0.031) — acetoin synthesis
- *entC*: 23, 19, 3 (q = 0.010) — siderophore biosynthesis
- No disease effect on any PGPR trait

### 4.3.4 Functional Redundancy

[Same as Paper 2, expanded with:]
- Redundancy distributions by compartment: Figure 4.4 (boxplots)
- Redundancy by health status: Figure 4.5
- Per-function redundancy variation: Figure 4.6 (identifying most and least redundant functions)

Key values:
- Bulk: 2.71 ± 1.60
- Endosphere: 2.99 ± 1.71
- Rhizosphere: 2.99 ± 1.67
- KW H = 122.59, p < 0.001
- Bulk vs endo: p < 0.001; bulk vs rhizo: p < 0.001; endo vs rhizo: p = 0.747
- Health: H = 0.037, p = 0.847 (NS)

### 4.3.5 Within-Group Functional Analysis

[Additional thesis-specific content:]
- Group-level redundancy: Table 4.6 (6 groups)
- Cross-compartment comparisons within same health state
- Confirmation that compartment signal dominates at all resolutions

---

## 4.4 Discussion

### 4.4.1 Compartment as the Master Variable for Functional Organisation

[Expanded from Paper 2 Discussion]
- >2,000 differential pathways vs 0 for disease
- Extends compositional finding (Chapter 3) to functional level
- Two-step selection model applies to function as well as taxonomy

### 4.4.2 Endosphere as a Nitrogen Fixation Niche

[Expanded from Paper 2]
- nifHDK enrichment matches maize (Zhang et al., 2022), rice (Edwards et al., 2015), poplar (Moyes et al., 2024)
- Aerobic N-fixation inside roots: physiological constraints and solutions
- Oil palm N demand and potential contribution of endophytic N-fixation
- Comparison of nif carrier counts with other systems

### 4.4.3 ACC Deaminase and Stress Mitigation

[Expanded from Paper 2]
- Endophyte benefit hypothesis: keeping host plant healthy benefits resident microbes
- acdS deletion impairs endophytic fitness (Nascimento et al., 2017)
- Relevance to BSR stress: ethylene as disease response hormone, acdS as moderator

### 4.4.4 Phosphate Cycling in Bulk Soil vs Root Compartments

[Expanded from Paper 2]
- gcd enrichment in bulk soil: mineral P solubilisation as free-living strategy
- Implications for P management in oil palm plantations
- Comparison with rhizosphere P-cycling in other crops

### 4.4.5 Functional Redundancy as Insurance

[Expanded from Paper 2]
- Root compartments maintain higher backup capacity
- Disease does NOT reduce redundancy — the insurance mechanism holds
- Comparison with redundancy patterns in other soil systems
- Ecological insurance hypothesis (Yachi & Loreau, 1999) supported
- Caution from Allison & Martiny (2008): redundancy is not universal

### 4.4.6 Limitations

[Same as Paper 2, expanded:]
- HQ-only annotation biases
- Gene presence ≠ gene expression
- Single timepoint
- Pathway completeness assumptions

---

## 4.5 Conclusions

[Mapped explicitly to Objective 2:]
1. Compartment drives massive functional differentiation (>2,000 pathways); disease has zero effect
2. Endosphere is enriched for N-fixation and stress mitigation (nifHDK, acdS)
3. Functional redundancy is higher in root compartments and maintained under disease
4. These results support the functional layer of the insurance hypothesis

---

## 4.6 References

[Superset of Paper 2 functional references + additional thesis-specific citations]
