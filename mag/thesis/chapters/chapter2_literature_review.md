# Chapter 2: Literature Review

**Estimated length:** 40–50 pages | **Writing priority:** 3rd (April 2026)

---

## 2.1 Oil Palm Agriculture and *Ganoderma* Basal Stem Rot

### 2.1.1 Oil Palm Biology and Cultivation

*Elaeis guineensis* Jacq. (Arecaceae) is a monocotyledonous perennial palm native to West Africa, now cultivated across the humid tropics. Key biological features:

- Lifespan: 25–30 years in commercial cultivation
- Root system: Adventitious primary roots from stem base, extensive secondary and tertiary roots
- Root compartments: Bulk soil → rhizosphere → rhizoplane → endosphere (distinct ecological niches)
- Annual yield: 3.5–4.0 t oil/ha (highest of any oilseed crop)

[EXPAND: Detailed root anatomy, compartment definitions with diagrams, rhizodeposition and exudate chemistry, palm physiology relevant to microbiome interactions]

### 2.1.2 *Ganoderma boninense* — Biology and Pathogenesis

- White-rot basidiomycete; lignin and cellulose degradation via laccase and peroxidase enzymes
- Infection primarily through root contact with infected stumps or soil inoculum
- Hemibiotrophic lifestyle: initial biotrophic colonisation followed by necrotrophic phase
- Fruiting bodies (basidiocarps) at basal stem indicate advanced infection

[EXPAND: Fungal biology, genetic diversity of Malaysian isolates, infection mechanism, virulence factors, host-pathogen molecular interactions]

### 2.1.3 BSR Epidemiology and Management

- Disease incidence increases across replanting cycles (inoculum buildup)
- Risk factors: soil type, previous cropping history, planting density
- Current management: cultural (sanitation, isolation), chemical (hexaconazole), biological (*Trichoderma*, AMF, *Pseudomonas*)
- Limitations of current approaches: No complete resistance, late detection, variable biocontrol efficacy

[EXPAND: Epidemiological studies across Malaysian estates, replanting cycle data, economic modelling, management comparison table]

---

## 2.2 Plant-Associated Microbiomes

### 2.2.1 The Root Microbiome Concept

- Definition: Microbial communities associated with root compartments (rhizosphere, rhizoplane, endosphere)
- Distinction from rhizosphere "effect" (Hiltner, 1904)
- Modern concept: holobiont framework (plant + microbiome as ecological unit)

[EXPAND: Historical development of rhizosphere concept, holobiont theory, plant-microbiome coevolution]

### 2.2.2 Compartment Concept and the Two-Step Selection Model

The two-step selection model (Bulgarelli et al., 2012; Edwards et al., 2015; Trivedi et al., 2020) proposes:

**Step 1 — Rhizosphere selection:** Root exudates (sugars, organic acids, amino acids) create a nutrient-rich zone that selects for a subset of bulk soil microorganisms. This is primarily a nutritional filter.

**Step 2 — Endosphere selection:** Plant immune responses (pattern-triggered immunity, effector-triggered immunity) and cell wall barriers further filter the community, admitting only organisms that can evade or suppress plant defences.

Evidence from model systems:
- *Arabidopsis*: Bulgarelli et al. (2012), Lundberg et al. (2012) — rhizosphere R² ~0.15–0.20, endosphere R² ~0.30–0.40
- Rice: Edwards et al. (2015) — compartment as largest source of variation
- Maize: Walters et al. (2018) — heritable microbiome components
- Barley: Bulgarelli et al. (2015) — wild vs domesticated comparison

[EXPAND: Detailed mechanism of each selection step, molecular basis of immune-microbe interactions, exudate chemistry, comparative table across plant species]

### 2.2.3 Tropical Crop Microbiomes

- Limited MAG-level studies compared with temperate model plants
- Existing work: sugarcane (de Souza et al., 2016), cassava, rubber
- Oil palm: 16S-based studies only (Shamsilawani et al., 2020)
- Gap: No genome-resolved analysis of tropical perennial root microbiomes

[EXPAND: Review of tropical crop microbiome studies, comparison of diversity patterns, effect of tropical soil conditions on microbiome assembly]

---

## 2.3 Metagenomics and MAGs

### 2.3.1 From Amplicons to Metagenomes

- 16S rRNA limitations: genus-level resolution, primer bias, no function
- Shotgun metagenomics: all DNA, functional annotation possible
- MAGs: genome-level resolution from shotgun data

[EXPAND: Technology comparison table, cost considerations, statistical power analysis]

### 2.3.2 MAG Recovery Pipeline

1. **Quality filtering** (fastp, Trimmomatic)
2. **Assembly** (MEGAHIT, metaSPAdes; co-assembly vs per-sample)
3. **Binning** (MetaBAT2, MaxBin2, CONCOCT, SemiBin; ensemble with DAS Tool)
4. **Quality assessment** (CheckM, CheckM2)
5. **Taxonomic classification** (GTDB-Tk)
6. **Abundance estimation** (CoverM, read mapping)

[EXPAND: Detailed pipeline flowchart, algorithm comparison, parameter selection rationale]

### 2.3.3 MIMAG Standards

Minimum Information about a Metagenome-Assembled Genome (Bowers et al., 2017):

| Quality tier | Completeness | Contamination | Additional requirements |
|---|---|---|---|
| High quality (HQ) | ≥90% | ≤5% | 23S, 16S, 5S rRNA; ≥18 tRNAs |
| Medium quality (MQ) | ≥50% | ≤10% | — |
| Low quality (LQ) | <50% | <10% | — |

[EXPAND: Quality assessment methodology, CheckM vs CheckM2, implications for downstream analysis]

---

## 2.4 Community Ecology of Soil Microbiomes

### 2.4.1 Alpha and Beta Diversity

- Shannon diversity, Simpson diversity, richness, evenness
- Bray-Curtis dissimilarity for compositional data
- Challenges of compositionality in relative abundance data

[EXPAND: Mathematical definitions, interpretation guidelines, compositional data analysis (CoDA) approaches]

### 2.4.2 Statistical Testing of Community Differences

- **PERMANOVA** (Anderson, 2001): Partitions total variation; R² quantifies effect size; pseudo-F tests significance
- **ANOSIM** (Clarke, 1993): Rank-based complementary test; R statistic interpretation
- **PERMDISP** (Anderson, 2006): Essential companion to PERMANOVA; tests dispersion homogeneity
- Pairwise comparisons with FDR correction

[EXPAND: Mathematical formulations, assumptions, power analysis, guidance on when each test is appropriate]

### 2.4.3 Compositional Data Analysis

- Compositional nature of relative abundance data (Aitchison, 1986; Gloor et al., 2017)
- Centre log-ratio (CLR) transformation
- Implications for differential abundance, correlation, and ordination

[EXPAND: CLR mathematics, comparison with other normalisation approaches, practical implementation]

### 2.4.4 Indicator Species and Core Microbiome

- IndVal method (Dufrêne & Legendre, 1997): specificity × fidelity
- Core microbiome: prevalence-based definition, threshold selection
- Compartment specificity: concentration index

[EXPAND: Statistical framework for IndVal, permutation testing, comparison with LEfSe and DESeq2 for biomarker identification]

---

## 2.5 Functional Annotation of MAGs

### 2.5.1 DRAM (Distilled and Refined Annotation of Metabolism)

- Multi-database integration: KEGG, UniRef90, Pfam, dbCAN, MEROPS (Shaffer et al., 2020)
- Advantages: standardised pipeline, reproducibility, community-level summaries
- KO (KEGG Orthology) as the functional unit

[EXPAND: DRAM algorithm, database versions, annotation confidence, comparison with other tools (Prokka, eggNOG)]

### 2.5.2 KEGG Pathway Analysis

- Hierarchical organisation: KO → pathway → module → category
- Pathway-level abundance: sum of carrier MAG abundances per sample
- Differential analysis: CLR transformation + non-parametric testing

[EXPAND: KEGG hierarchy, pathway completeness estimation, interpretation caveats]

### 2.5.3 CAZyme Classification

- Six CAZyme classes: GH, GT, CE, AA, CBM, PL (Lombard et al., 2014)
- dbCAN for family-level classification
- Ecological relevance: carbon cycling, cell wall degradation, exopolysaccharide production

[EXPAND: CAZyme family descriptions with ecological roles, relevance to soil organic matter cycling, plant cell wall interactions]

### 2.5.4 Plant Growth-Promoting Traits

Key PGPR traits and their molecular markers:

| Trait | Gene(s) | Mechanism | Ecological role |
|-------|---------|-----------|-----------------|
| N fixation | *nifH/D/K* | Nitrogenase complex | Direct N supply to plant |
| P solubilisation | *gcd*, *pqqC* | Glucose dehydrogenase, PQQ cofactor | Mineral P release |
| ACC deaminase | *acdS* | Ethylene precursor cleavage | Stress reduction |
| IAA biosynthesis | *ipdC*, *acsA* | Indole-3-acetic acid production | Root growth promotion |
| Siderophore | *entA/B/C* | Iron chelation | Fe acquisition, pathogen competition |
| Acetoin/butanediol | *budB* | Volatile organic compound | ISR induction |
| Phenazine | *phzF* | Antibiotic production | Direct pathogen inhibition |

[EXPAND: Molecular mechanisms for each trait, literature on endosphere vs rhizosphere enrichment patterns, relevance to BSR suppression]

### 2.5.5 Functional Redundancy

- Definition: Multiple taxa carrying the same function (insurance against species loss)
- Quantification: Shannon diversity of function carriers
- Ecological insurance hypothesis (Yachi & Loreau, 1999)
- Evidence for and against redundancy in soil systems (Louca et al., 2016; Allison & Martiny, 2008)

[EXPAND: Mathematical framework, redundancy vs complementarity, relationship to community stability]

---

## 2.6 Co-occurrence Network Analysis

### 2.6.1 Rationale for Network Approaches

- Pairwise correlations capture co-occurrence patterns (potential interactions)
- Networks reveal emergent properties: modularity, keystones, hubs/bridges
- Compositional data requires specialised correlation measures

[EXPAND: Network theory primer, distinction between co-occurrence and interaction, ecological interpretation guidelines]

### 2.6.2 Phi Proportionality

- Standard correlation (Pearson, Spearman) on relative abundance produces spurious correlations (Pearson, 1897; Gloor et al., 2017)
- Phi proportionality (Lovell et al., 2015): measures proportionality between components in compositional data
- Phi is invariant to total library size, avoiding false positive associations

[EXPAND: Mathematical derivation, comparison with SparCC, SPIEC-EASI, relationship to Aitchison log-ratio distance]

### 2.6.3 Modularity and Community Detection

- Modularity (*Q*; Newman, 2006): quantifies how well a network partitions into modules
- Louvain algorithm (Blondel et al., 2008): fast, scalable community detection
- Null model validation: observed *Q* vs random expectation (z-score)
- Threshold sensitivity: testing network robustness across edge-filtering stringencies

[EXPAND: Modularity mathematics, algorithm comparison, null model construction, interpretation of z-scores]

### 2.6.4 Keystone Taxa

- Definition: Organisms whose removal disproportionately alters community structure (Power et al., 1996)
- Network-based identification: high centrality (betweenness, closeness, degree) + low abundance
- Composite scoring: normalised centrality inverse-weighted by abundance

[EXPAND: Keystone concept history, detection methods comparison, validation approaches, caveats]

### 2.6.5 Hub-Bridge Classification (Guimerà-Amaral Framework)

- Within-module degree (*z*): connectivity within own module
- Participation coefficient (*P*): evenness of connections across modules
- Four roles: hub (*z* > 2.5), bridge (*P* > 0.62), peripheral, kinless
- Ecological interpretation: hubs maintain local stability, bridges maintain global connectivity

[EXPAND: Mathematical definitions, ecological role of each classification, examples from soil microbiome literature]

### 2.6.6 Differential Network Analysis

- Comparing networks between conditions (compartments, health states)
- Edge classification: conserved, gained, lost
- Interpretation: rewiring vs densification vs fragmentation

[EXPAND: Methods for network comparison, statistical testing of differential edges, biological interpretation]

---

## 2.7 Microbial Ecology of Oil Palm Soils

### 2.7.1 Previous 16S rRNA Studies

- Shamsilawani et al. (2020): Rhizosphere and endosphere of healthy vs *Ganoderma*-infected oil palm
  - Found bacterial community shifts associated with disease
  - Limited to genus-level resolution; no functional or network analysis
- Other oil palm soil studies: limited to bulk soil bacterial/fungal surveys

[EXPAND: Comprehensive review of all published oil palm microbiome studies, comparison table, identified gaps]

### 2.7.2 Tropical Crop Disease Microbiomes

- Fusarium wilt of banana: network complexity changes (Deng et al., 2021; Gao et al., 2021)
- Verticillium wilt: Myxococcota enrichment (Zhao et al., 2024)
- Rice blast: endosphere microbiome responses
- Common patterns: disease-induced structural changes, variable compositional effects

[EXPAND: Cross-system comparison of disease effects on microbiome structure, universal vs pathosystem-specific patterns]

---

## 2.8 Myxococcota and Acidobacteriota Ecology

### 2.8.1 Myxococcota as Soil Micropredators

- >60% of soil bacterivores are Myxococcota (Petters et al., 2021)
- Predation mechanism: swarming, lytic enzyme secretion, prey cell lysis, nutrient uptake
- Dominant predatory families: Haliangiaceae, Polyangiaceae (class Polyangia)
- Prey range: Gram-negative > Gram-positive (Morgan et al., 2010)
- Biocontrol potential: *Corallococcus* sp. EGB reduced Fusarium wilt 54–80% (Ye et al., 2020)

[EXPAND: Myxococcota lifecycle, predation mechanisms, secretion systems (T3SS/T6SS), secondary metabolite diversity, biocontrol literature, phylogenetic diversity in soil]

### 2.8.2 Acidobacteriota as Keystone Decomposers

- Ubiquitous in soil: 5–70% of communities (Kielak et al., 2016)
- Predominantly oligotrophic: slow growth, recalcitrant carbon specialists
- *Terracidiphilus gabretensis*: degrades cellulose, chitin, CMC; 132 active genes in soil (García-Fraile et al., 2016)
- Keystone function: carbon cycling, exopolysaccharide production, nitrite utilisation (Kalam et al., 2020)
- Plant growth-promoting potential: N-fixation, P-solubilisation, siderophores (Goncalves et al., 2024)

[EXPAND: Acidobacteriota subdivisions, ecological strategies (oligotrophy vs copiotrophy), CAZyme repertoire, role in soil organic matter cycling]

### 2.8.3 The Acidobacteriota-Myxococcota Trophic Link

- Negative correlation between abundances (Zhou et al., 2020: r = −0.495, p < 0.001)
- Proposed predator-prey relationship
- Nutrient recycling: Acidobacteriota degrade carbon → Myxococcota prey on Acidobacteriota → release N, P, C
- Both associated with nutrient cycling indices in farmland soils (Dai et al., 2023)

[EXPAND: Trophic cascade model, evidence from multiple studies, comparison with marine microbial loops, implications for soil nutrient cycling]

---

## 2.9 Summary and Research Justification

This literature review identifies the following knowledge gaps that this thesis addresses:

1. **No MAG-level study of oil palm root microbiomes** exists, despite the economic importance of the crop and the severity of BSR
2. **Compartment-driven assembly** has not been tested in oil palm at genome resolution
3. **Functional capacity** (PGPR traits, functional redundancy) is unknown for oil palm root-associated MAGs
4. **Co-occurrence networks** and keystone taxa have not been characterised in oil palm soils
5. **The insurance hypothesis** (taxonomic generalism + functional redundancy + modular networks = resilience) has not been tested in any crop disease system using integrated MAG-level analyses

These gaps motivate the three research objectives outlined in Chapter 1.

---

## References

[To be compiled — estimated 200–300 references covering all sections above]
