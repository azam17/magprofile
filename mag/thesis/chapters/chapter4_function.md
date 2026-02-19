# Chapter 4: Functional Capacity and Plant Growth-Promoting Potential of Compartment-Stratified MAGs

*This chapter addresses Objective 2 and corresponds to the functional profiling component of Paper 2.*

---

## 4.1 Introduction

[Expanded from Paper 2 Introduction, focused specifically on functional questions. ~3-4 pages]

Key points:
- Chapter 3 findings raise functional questions that taxonomy alone cannot answer
- Table linking each ecological finding to a specific functional question (from FINDINGS_v021.md Section 3):
  - Myxococcota guild → lytic enzymes, secretion systems?
  - Roseiarcus in roots → photosynthetic genes retained?
  - Acidobacteriota rhizosphere → CAZyme enrichment?
  - Archaea exclusion → nitrogen cycling restructured?
  - Verrucomicrobiota disease markers → metabolic traits?
  - Ultrasmall parasites → host dependency genes?
- Functional redundancy as an insurance mechanism
- Specific hypotheses for this chapter

## 4.2 Materials and Methods

### 4.2.1 MAG selection for functional annotation

- 567 HQ MAGs (completeness ≥90%, contamination ≤5%) selected from 1,018 total
- Rationale for excluding MQ MAGs: gene absence on incomplete genomes is unreliable (a MAG at 60% completeness has 40% chance of missing any gene)
- Statistical power trade-off acknowledged

### 4.2.2 DRAM annotation pipeline

**Source:** Paper 2 Methods, v043_analysis_report.txt Section 2.1

- DRAM v1.x pipeline configuration
- Database versions (KEGG, UniRef90, Pfam, dbCAN, MEROPS)
- Gene calling parameters (Prodigal)
- KO assignment methodology
- Quality filtering and confidence thresholds

### 4.2.3 Functional profiling analyses

**Pathway abundance calculation:**
- Per-sample pathway abundance = sum of relative abundances of MAGs carrying each KO
- CLR transformation for differential testing

**CAZyme classification:**
- dbCAN annotation into 6 classes (GH, GT, CE, AA, CBM, PL)
- Family-level and class-level abundance quantification

**PGPR trait screening:**
- 13-gene panel: nifH, nifD, nifK, gcd, pqqC, acdS, ipdC, acsA, budB, entA, entB, entC, phzF
- Presence/absence determination from KO annotations
- Trait prevalence comparison across compartments (Kruskal-Wallis, FDR)

**Functional redundancy:**
- Shannon diversity of MAGs carrying each function, per sample
- Comparison across compartments (KW) and health states (Mann-Whitney)

**Differential functional analysis:**
- Mann-Whitney U on CLR-transformed abundances
- FDR correction at q < 0.05
- All pairwise comparisons: compartment (3), health (1), group (15)

## 4.3 Results

### 4.3.1 Annotation overview

**Source:** v043_analysis_report.txt Section 2.1

- 6,233 KEGG orthologues (KOs) detected
- 532 CAZyme families across 6 classes
- 13 PGPR traits screened
- 567 MAGs annotated

CAZyme distribution (Table 4.1):

| Class | Genes | Families | MAGs |
|-------|------:|--------:|-----:|
| GH | 26,069 | 193 | 557 (98.2%) |
| GT | 22,083 | 97 | 567 (100%) |
| CE | 5,155 | 42 | 529 (93.3%) |
| AA | 3,764 | 10 | 530 (93.5%) |
| CBM | 2,061 | 159 | 438 (77.2%) |
| PL | 824 | 31 | 301 (53.1%) |

[Thesis addition: KO distribution histogram, MAG-level annotation completeness, comparison with published soil metagenome annotations]

### 4.3.2 Massive compartment-driven functional differentiation

**Source:** v043_analysis_report.txt Section 2.4

- Bulk vs Endosphere: ~2,000+ significant KO pathways (q < 0.05)
- Bulk vs Rhizosphere: ~1,500+ significant pathways
- Endosphere vs Rhizosphere: ~500+ significant pathways
- CAZyme differentials: ~50-200+ per comparison

Disease comparisons:
- Diseased vs Healthy: 0 significant pathways (all q > 0.99)
- 0 significant CAZymes (all q > 0.47)
- Within-compartment disease comparisons: all 0

[Thesis addition: Volcano plots for each comparison, top enriched/depleted pathways by compartment, KEGG module completeness analysis]

### 4.3.3 PGPR trait enrichment by compartment

**Source:** v043_analysis_report.txt Sections 2.3, 5.1–5.4

**Nitrogen fixation (nifHDK):**
- nifH: 0 bulk, 5 endo, 1 rhizo (q = 0.031)
- nifD: 0, 4, 1 (p = 0.042, q = 0.078 — trend)
- nifK: 0, 6, 1 (q = 0.016)
- Coherent biological signal: all three subunits enriched in endosphere

**ACC deaminase (acdS):**
- 5 bulk, 15 endo, 3 rhizo (q = 0.006)
- Endophytic stress relief mechanism

**Phosphate solubilisation (gcd):**
- 120 bulk, 117 endo, 84 rhizo (q = 0.006)
- Free-living P-cycling dominates

**Other significant traits:**
- budB (acetoin): 162, 131, 103 (q = 0.031)
- entC (siderophore): 23, 19, 3 (q = 0.010)

**No health status effects on any PGPR trait.**

[Thesis addition: Carrier MAG taxonomy for each trait, cross-tabulation with Chapter 3 indicator species, comparison with maize/rice PGPR literature]

### 4.3.4 Functional redundancy: higher in roots, unaffected by disease

**Source:** v043_analysis_report.txt Section 2.2

By compartment:
- Bulk: 2.71 ± 1.60
- Endosphere: 2.99 ± 1.71
- Rhizosphere: 2.99 ± 1.67
- Kruskal-Wallis H = 122.59, p < 0.001
- Mann-Whitney bulk vs endo: p < 0.001
- Mann-Whitney bulk vs rhizo: p < 0.001
- Mann-Whitney endo vs rhizo: p = 0.747 (NS)

By health status:
- Diseased: 3.00 ± 1.70
- Healthy: 2.99 ± 1.69
- KW H = 0.037, p = 0.847 (NS)

By group:
- KW H = 264.03, p < 0.001
- Pattern driven entirely by compartment component

[Thesis addition: Redundancy distribution histograms, function-level redundancy breakdown, comparison with Louca et al. (2016) and Allison & Martiny (2008)]

### 4.3.5 MAG count discrepancy (567 vs 1,018)

**Source:** v043_analysis_report.txt Section 5.4

- 451 MQ MAGs excluded from functional annotation
- Rationale: 60% complete MAG has 40% chance of missing any given gene
- Trade-off: reduced power for rare traits (nif genes: 5-7 carriers out of 567)
- Defensible: gene absence on incomplete genomes is unreliable

## 4.4 Discussion

### 4.4.1 Compartment controls function as well as composition

[Expanded from Paper 2 Discussion. ~3 pages]

- Two-step selection model extends to function
- >2,000 differential pathways = comprehensive functional reshaping
- Zero disease effect on any functional metric — remarkable resilience
- Comparison with Edwards et al. (2015) rice, Zhang et al. (2022) maize

### 4.4.2 Endosphere as a nitrogen fixation niche

[Expanded from Paper 2. ~2 pages]

- nifHDK coherent enrichment despite low carrier count
- Literature: Zhang et al. (2022) nifH in maize xylem, Moyes et al. (2024) poplar endophyte N-fixation
- Low-oxygen endosphere environment may favour nitrogenase activity
- Limitation: gene presence ≠ activity, metatranscriptomic validation needed

### 4.4.3 ACC deaminase: endophytes reduce host stress

[Expanded from Paper 2. ~1 page]

- acdS as endosphere colonisation fitness factor (Nascimento et al., 2017)
- Glick (2014): ACC deaminase as key PGPR mechanism
- Relevance to BSR stress response

### 4.4.4 Functional redundancy as insurance

[Expanded from Paper 2. ~2 pages]

- Insurance hypothesis (Yachi and Loreau, 1999)
- Higher redundancy in root compartments = more backup capacity
- Disease does not erode redundancy = insurance mechanism is robust
- Hemkemeyer et al. (2024): plant-driven rhizosphere effects maintain redundancy
- Allison & Martiny (2008): redundancy is not universal — this finding is notable

### 4.4.5 Implications for BSR management

[~1 page]

- PGPR-carrying endophytes as biocontrol candidates
- Functional resilience suggests microbiome augmentation rather than replacement
- Targeted inoculants should carry nifHDK + acdS for endosphere colonisation

## 4.5 Summary

- Compartment drives massive functional differentiation (>2,000 DE pathways)
- Disease has zero effect on function, PGPR traits, or redundancy
- Endosphere is enriched for nitrogen fixation and ACC deaminase
- Root compartments maintain higher functional redundancy than bulk soil
- The functional insurance layer is preserved under BSR

*These functional findings complement the ecological patterns of Chapter 3 and motivate the network analysis in Chapter 5: if both composition and function are compartment-driven and disease-resilient, does the interaction architecture show the same pattern?*

---

*[Data sources: v043_analysis_report.txt Sections 2.1–2.4, 5.1–5.4]*
