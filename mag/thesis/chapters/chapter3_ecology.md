# Chapter 3: Compartment-Driven Assembly of Oil Palm Root Microbiomes under Ganoderma BSR

**Estimated length:** 30–40 pages | **Writing priority:** 4th (May 2026)

**Relationship to Paper 1:** This chapter expands Paper 1 (FEMS Microbiology Ecology) into thesis format. Additional content includes: expanded methods, supplementary results, extended discussion with more literature, and thesis-specific framing connecting to Objectives 1–3.

**Source:** Paper 1 manuscript (`thesis/papers/paper1_ecology.md`) + `FINDINGS_v021.md` + `v043_analysis_report.txt` Section 10

---

## 3.1 Introduction

[Adapted from Paper 1 Introduction, expanded with:]
- Explicit link to Objective 1
- Extended context on compartment ecology in perennial crops
- Deeper coverage of BSR-microbiome interactions literature
- Justification for MAG-level approach over 16S

---

## 3.2 Materials and Methods

### 3.2.1 Study Site and Experimental Design

[Same as Paper 1, with additional detail:]
- Plantation coordinates, soil type, climate data
- Palm age, variety, management history
- BSR confirmation criteria (fruiting body + DIS score)
- Sampling schematic diagram (Figure 3.1)

### 3.2.2 DNA Extraction and Metagenomic Sequencing

[Same as Paper 1, with additional detail:]
- DNA quality metrics (Nanodrop, Qubit, gel electrophoresis)
- Library size distributions
- Sequencing depth per sample (total reads, after filtering)

### 3.2.3 Metagenome Assembly and MAG Recovery

[Same as Paper 1, with additional detail:]
- Co-assembly strategy rationale
- Per-algorithm bin counts before DAS Tool refinement
- Dereplication parameters (ANI threshold)
- CheckM2 quality distribution plots

### 3.2.4 Taxonomic Classification and Abundance Estimation

[Same as Paper 1]

### 3.2.5 Ecological Analyses

[Same as Paper 1, with additional detail on:]
- Alpha diversity: justification of metrics, rarefaction considerations
- Beta diversity: distance metric selection, PERMANOVA model specification
- CLR transformation: zero handling, pseudocount strategy
- Indicator species: permutation parameters, IndVal threshold selection
- Compartment specificity: index derivation, specialist threshold justification

### 3.2.6 Software and Reproducibility

[Same as Paper 1, plus:]
- Full command-line invocations
- Package version pinning
- Data availability statement

---

## 3.3 Results

### 3.3.1 MAG Recovery and Quality

- 1,018 non-redundant MAGs (567 HQ, 451 MQ)
- 36 phyla, 85 classes, 320 genera
- Quality distribution: Figure 3.2 (completeness vs contamination scatter)
- Taxonomic composition: Figure 3.3 (phylum-level barplot across samples)

### 3.3.2 Alpha Diversity

[Same as Paper 1 Results, Table 3.1]
- Shannon 5.99 ± 0.23, Simpson 0.994, richness 1,018
- No significant differences by compartment or health status
- Evenness range 0.787–0.921

### 3.3.3 Beta Diversity and Community Composition

[Same as Paper 1, expanded with:]
- PCoA ordination plots: Figure 3.4 (compartment), Figure 3.5 (health status)
- Complete PERMANOVA table: Table 3.2
- All pairwise PERMANOVA results: Table 3.3
- PERMDISP within-group distances: Table 3.4
- Variance partitioning: compartment 42.4% vs health 5.2% vs residual

### 3.3.4 Differential Abundance

[Same as Paper 1, expanded with:]
- Volcano plots for all 6 compartment pairwise comparisons: Figure 3.6
- Full count table: Table 3.5 (634/474/376 by compartment, 0 by health)
- Effect size distributions for top MAGs
- Cross-compartment comparisons within health states

### 3.3.5 Indicator Species

[Same as Paper 1, expanded with:]
- Extended top-20 indicator tables for each compartment: Tables 3.6–3.8
- Phylum composition of indicators per compartment: Figure 3.7
- Disease vs healthy indicator comparison: Table 3.9, Figure 3.8
- Complete disease indicator list (125 MAGs): Supplementary Table S1

### 3.3.6 Core Microbiome

[Same as Paper 1]
- Prevalence curve: Figure 3.9
- 1,017/1,018 at 100% prevalence across all groupings

### 3.3.7 Compartment Specificity

[Same as Paper 1, expanded with:]
- Specificity score distribution: Figure 3.10
- 65 endosphere specialists characterised by phylum: Table 3.10
- Generalist vs specialist comparison: abundance, prevalence, taxonomy

### 3.3.8 Endosphere Guild Discovery

#### 3.3.8.1 Myxococcota Predatory Guild

[Expanded from Paper 1]
- 15 MAGs, all Polyangia
- Genera: *Palsa-1150* (6), JADGRB01 (7), others (2)
- IndVal scores: 0.87–0.91
- 10-fold enrichment (4.6% endo vs 0.4% bulk, d = -3.86)
- Compartment abundance profiles per MAG: Figure 3.11

#### 3.3.8.2 *Roseiarcus* Photoheterotrophs

[Expanded from Paper 1]
- 5 MAGs among top 20 endosphere indicators
- IndVal 0.87–0.90
- Photoheterotrophic paradox discussion

#### 3.3.8.3 Archaea Exclusion

[Expanded from Paper 1]
- Thermoplasmatota: d = 6.43
- Thermoproteota: d = 6.02
- Barplot of archaeal phyla across compartments: Figure 3.12

### 3.3.9 Disease-Associated Taxa

[Expanded from Paper 1]
- 3 FDR-significant MAGs (2 Verrucomicrobiota, 1 Actinomycetota)
- Ultrasmall dysbiosis markers (Micrarchaeota, Patescibacteria)
- 125 vs 36 indicator asymmetry explained

---

## 3.4 Discussion

[Expanded from Paper 1 Discussion, with:]
- Extended comparison with other crop systems (rice, maize, barley, *Arabidopsis*, sugarcane)
- Deeper ecological theory: niche filtering, environmental filtering vs competitive exclusion
- Extended Myxococcota discussion: lifecycle, predation mechanisms, biocontrol potential
- *Roseiarcus* paradox: metabolic switching hypotheses
- Acidobacteriota rhizosphere dominance: oligotrophy and root exudate specialisation
- Archaea exclusion: comparison with other plant systems, mechanistic hypotheses
- Disease response: precision pathobiome concept, comparison with other disease systems
- Universal membership: implications for microbiome manipulation strategies

### 3.4.1 Compartment-Driven Assembly in Tropical Perennial Crops
### 3.4.2 Novel Endosphere Guilds: The Myxococcota Predatory Niche
### 3.4.3 The Roseiarcus Photoheterotrophic Paradox
### 3.4.4 Compositional Boundary at the Rhizoplane
### 3.4.5 Archaea Exclusion from Root Compartments
### 3.4.6 Fidelity-Based Disease Signals and Dysbiosis Markers
### 3.4.7 Universal Membership and the "How Much?" Paradigm
### 3.4.8 Limitations

---

## 3.5 Conclusions

[Same as Paper 1, with explicit mapping to Objective 1]

---

## 3.6 References

[Superset of Paper 1 references + additional thesis-specific citations]
