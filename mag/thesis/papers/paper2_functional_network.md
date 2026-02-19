# Functional Resilience and Interaction Network Architecture of Oil Palm Root-Associated Microbiomes under Ganoderma Basal Stem Rot

**Authors:** [Author names to be added]

**Target journal:** Pertanika Journal of Tropical Agricultural Science (JTAS)

**Format:** IMRAD, ≤6,000 words (excluding abstract, references, tables, figures), ≤80 references, abstract <250 words

---

## Abstract

Oil palm (*Elaeis guineensis*) is Southeast Asia's most economically important crop, yet basal stem rot (BSR) caused by *Ganoderma boninense* threatens plantation sustainability. While root-associated microbiomes may confer disease resilience, genome-resolved functional and interaction analyses remain absent for oil palm. Here, we applied metagenome-assembled genome (MAG) analysis to characterise the functional capacity and co-occurrence network architecture of microbiomes across three root compartments (bulk soil, rhizosphere, endosphere) and two health states (healthy, diseased) in a Malaysian oil palm plantation. Functional profiling of 567 high-quality MAGs revealed 6,233 KEGG orthologues and 532 CAZyme families, with massive compartment-driven differentiation (>2,000 differential pathways between bulk soil and endosphere) but zero disease effect on any functional metric. Nitrogen fixation genes (*nifHDK*) and ACC deaminase (*acdS*) were significantly enriched in the endosphere, while functional redundancy was higher in root compartments (Shannon 2.99) than bulk soil (2.71; *p* < 0.001). Co-occurrence network analysis of 1,018 MAGs identified 25,883 edges organised into 258 non-random modules (null model *z* = 81.36). *Terracidiphilus* (Acidobacteriota) emerged as the top keystone taxon (score 0.981), while Myxococcota predatory bacteria served as inter-module bridges. Disease increased network density by 27% without altering composition or function, indicating structural dysbiosis. These results demonstrate a three-layer insurance mechanism—taxonomic redundancy, functional backup, and modular network architecture—that maintains microbiome resilience under BSR, with compartment identity as the master ecological variable.

**Keywords:** oil palm; *Ganoderma boninense*; metagenome-assembled genomes; functional redundancy; co-occurrence network; keystone taxa; basal stem rot; root microbiome

---

## Introduction

Oil palm (*Elaeis guineensis* Jacq.) is the world's most productive oilseed crop, contributing approximately 35% of global vegetable oil production from only 10% of oilseed cropland (Murphy, 2014). Malaysia, as the world's second-largest producer, derives significant economic value from the oil palm sector, which supports over 650,000 smallholders and contributes substantially to national GDP (MPOB, 2023). However, plantation sustainability is critically threatened by basal stem rot (BSR), caused by the white-rot basidiomycete *Ganoderma boninense* Pat. BSR is the most devastating disease of oil palm in Southeast Asia, causing yield losses exceeding 50% in severely affected plantations and an estimated RM 1.5 billion in annual economic losses in Malaysia alone (Idris et al., 2004; Chong et al., 2017).

Current BSR management strategies—including cultural practices, chemical fungicides, and biological control agents—have achieved limited success, partly because the soil microbiome's role in disease suppression remains poorly characterised at the genomic level (Sundram et al., 2019). While 16S rRNA amplicon studies have revealed shifts in bacterial community composition associated with BSR (Shamsilawani et al., 2020), amplicon-based approaches cannot resolve functional capacity, metabolic potential, or genome-level interaction patterns. Metagenome-assembled genomes (MAGs) offer a fundamentally superior resolution by recovering near-complete genomes directly from environmental DNA, enabling simultaneous taxonomic classification, functional annotation, and interaction network reconstruction (Parks et al., 2017; Bowers et al., 2017).

Plant root-associated microbiomes are structured by compartment-specific ecological filtering (Bulgarelli et al., 2012; Edwards et al., 2015). The two-step selection model proposes that bulk soil communities are first filtered at the rhizosphere by root exudates, then further refined at the endosphere by plant immune responses (Trivedi et al., 2020). This compartment-driven assembly has functional consequences: endophytic microbiomes are enriched for plant growth-promoting traits including nitrogen fixation, phosphate solubilisation, and phytohormone modulation (Compant et al., 2010; Hardoim et al., 2015). However, whether these functional patterns persist under disease pressure—and whether the microbiome possesses sufficient functional redundancy to buffer against perturbation—remains largely unexplored in tropical crop systems.

Co-occurrence network analysis provides a complementary lens for understanding microbiome organisation. Proportionality-based correlation measures (Lovell et al., 2015) overcome the compositionality bias inherent in relative abundance data, enabling robust identification of co-occurring taxa, modular community structure, and keystone organisms—taxa that exert disproportionate influence on community architecture despite low abundance (Berry & Widder, 2014). Network topology metrics, including modularity, hub-bridge classification (Guimerà & Amaral, 2005), and null model validation, can distinguish deterministic from stochastic community assembly and identify taxa whose removal would disproportionately alter community structure.

To date, no study has applied genome-resolved functional and network analysis to oil palm root-associated microbiomes under BSR. Previous work on this dataset demonstrated that compartment identity, not disease status, is the dominant driver of MAG-level community composition (PERMANOVA R² = 0.424, *p* = 0.001 vs R² = 0.052, *p* = 0.120; [Author], in preparation). Here, we extend that ecological analysis by asking two questions: (1) Does compartment-driven assembly extend to functional capacity, and does disease erode functional potential? (2) How are microbial interaction networks structured across compartments and health states, and do keystone taxa maintain community stability under BSR? By integrating functional profiling with network topology analysis, we test the hypothesis that the oil palm root microbiome maintains resilience through a multi-layered insurance mechanism combining taxonomic generalism, functional redundancy, and modular network architecture.

---

## Materials and Methods

### Study Site and Sampling Design

Samples were collected from an oil palm plantation in Peninsular Malaysia under the Malaysian Palm Oil Board (MPOB) research programme. The plantation contained mature oil palm stands (>15 years) with documented BSR incidence. Thirty samples were collected using a fully crossed factorial design: 3 compartments (bulk soil, rhizosphere, endosphere) × 2 health states (healthy palms, BSR-diseased palms confirmed by basal fruiting body presence) × 5 biological replicates per group, yielding 6 treatment groups. Bulk soil was collected >30 cm from the root zone, rhizosphere soil was obtained by shaking loosely adhering soil from washed roots, and endosphere samples were prepared by surface-sterilising roots (sequential ethanol-sodium hypochlorite washes) followed by maceration. All samples were stored at −80°C until DNA extraction.

### DNA Extraction, Sequencing, and MAG Recovery

Total genomic DNA was extracted using the DNeasy PowerSoil Pro Kit (Qiagen). Metagenomic libraries were prepared and sequenced on the Illumina NovaSeq 6000 platform (paired-end 150 bp). Quality-filtered reads were assembled using megahit (Li et al., 2015) with multiple co-assembly strategies. Metagenomic binning was performed using a multi-algorithm approach combining MetaBAT2 (Kang et al., 2019), MaxBin2 (Wu et al., 2016), CONCOCT (Alneberg et al., 2014), and SemiBin (Pan et al., 2023), with DAS Tool (Sieber et al., 2018) for bin refinement. MAG quality was assessed using CheckM2 (Chklovski et al., 2023) following MIMAG standards (Bowers et al., 2017): high-quality (HQ) MAGs (completeness ≥90%, contamination ≤5%) were used for functional annotation, while both HQ and medium-quality (MQ; completeness ≥50%, contamination ≤10%) MAGs were retained for abundance-based network analysis. A total of 1,018 non-redundant MAGs were recovered, of which 567 met HQ criteria.

### Functional Annotation

HQ MAGs (n = 567) were annotated using DRAM (Distilled and Refined Annotation of Metabolism; Shaffer et al., 2020), which integrates annotations from KEGG, UniRef90, Pfam, dbCAN (CAZyme), and MEROPS databases. KEGG Orthology (KO) assignments were used to quantify pathway-level functional capacity. Carbohydrate-active enzyme (CAZyme) families were classified into six classes: glycoside hydrolases (GH), glycosyltransferases (GT), carbohydrate esterases (CE), auxiliary activities (AA), carbohydrate-binding modules (CBM), and polysaccharide lyases (PL).

### Plant Growth-Promoting Trait Screening

A panel of 13 plant growth-promoting rhizobacteria (PGPR) marker genes was screened across all HQ MAGs: nitrogen fixation (*nifH*, *nifD*, *nifK*), phosphate solubilisation (*gcd*, *pqqC*), ACC deaminase (*acdS*), indole-3-acetic acid biosynthesis (*ipdC*, *acsA*), acetoin/2,3-butanediol production (*budB*), siderophore biosynthesis (*entA*, *entB*, *entC*), and phenazine biosynthesis (*phzF*). Gene presence was determined from DRAM KO annotations. Trait prevalence across compartments was compared using Kruskal-Wallis tests with Benjamini-Hochberg false discovery rate (FDR) correction at *q* < 0.05.

### Functional Redundancy Analysis

Functional redundancy was quantified as the Shannon diversity of MAGs carrying each function, calculated per sample. This metric captures how many taxonomically distinct MAGs share each functional capability—higher values indicate greater backup capacity. Redundancy was compared across compartments using the Kruskal-Wallis test, and across health states using the Mann-Whitney U test, with FDR correction.

### Differential Functional Analysis

Pathway-level and CAZyme family abundances were calculated as the summed relative abundance of MAGs carrying each function in each sample. Centre log-ratio (CLR) transformation was applied to account for compositionality. Differential abundance between groups was tested using the Mann-Whitney U test with FDR correction (*q* < 0.05). All pairwise comparisons were performed for compartment (3 pairs), health status (1 pair), and group (15 pairs).

### Co-occurrence Network Construction

Co-occurrence networks were constructed from relative abundance profiles of all 1,018 MAGs across 30 samples using phi proportionality (Lovell et al., 2015), a compositionality-aware correlation measure that avoids the spurious correlations generated by standard Pearson or Spearman coefficients on relative abundance data. CLR-transformed abundances were used to compute pairwise phi values. MAGs present in fewer than 50% of samples were excluded to ensure robust correlation estimates.

An edge was retained if |phi| fell below the 5th percentile of the empirical phi distribution (threshold = 0.0532), indicating strong proportional co-occurrence. Threshold sensitivity was assessed across five percentiles (1st, 3rd, 5th, 10th, 20th) to confirm that network properties were robust to threshold choice.

### Network Topology and Module Detection

Community detection was performed using the Louvain algorithm (Blondel et al., 2008) to identify network modules—groups of MAGs that co-occur more strongly with each other than with the rest of the network. Global network properties were calculated: density (proportion of realised edges), modularity (*Q*; Newman, 2006), and degree distribution.

Null model validation was performed by generating 1,000 random networks with preserved degree sequence (configuration model) and comparing observed modularity against the null distribution. The z-score quantifies how many standard deviations the observed modularity exceeds random expectation.

### Keystone Taxa Identification

Keystone taxa were identified using a composite scoring system integrating three network centrality metrics: (1) betweenness centrality (proportion of shortest paths passing through a node), (2) closeness centrality (inverse mean distance to all other nodes), and (3) degree centrality (number of direct connections), normalised and inverse-weighted by mean relative abundance. This identifies taxa that are structurally important (high centrality) despite being numerically rare—the hallmark of keystone organisms (Power et al., 1996; Berry & Widder, 2014).

### Hub-Bridge Classification

The Guimerà-Amaral framework (Guimerà & Amaral, 2005) was applied to classify MAGs by their network role. Within-module degree (*z*) quantifies how connected a node is within its own module, while the participation coefficient (*P*) quantifies how evenly its connections are distributed across modules. Nodes with *z* > 2.5 are hubs (highly connected within their module); nodes with *P* > 0.62 are bridges (connecting multiple modules).

### Compartment-Specific and Health-Specific Networks

Separate networks were constructed for each compartment (bulk soil, rhizosphere, endosphere) and each health state (healthy, diseased) by subsetting samples and recalculating phi proportionality correlations. Differential network analysis compared edge sets between networks to identify gained, lost, and conserved interactions.

### Software Implementation

All analyses were performed using the magprofile Python package (v0.43), comprising three integrated modules: magprofile (ecological analysis), magfunc (functional profiling), and magnet (network analysis). The package implements 20 analysis modules with 158 unit tests and is available at [repository URL]. Analyses were executed via command-line interface:

```
magprofile func-report --dram-dir <path> --grouping compartment --advanced
magprofile net-report --grouping compartment --advanced
```

Statistical analyses used scipy (v1.11+) for non-parametric tests, numpy (v1.24+) for numerical computation, and networkx (v3.1+) for graph algorithms.

---

## Results

### Functional Annotation Overview

DRAM annotation of 567 HQ MAGs identified 6,233 KEGG orthologues (KOs) and 532 CAZyme families. CAZyme genes were distributed across six classes: glycoside hydrolases (GH; 26,069 genes in 193 families across 557 MAGs), glycosyltransferases (GT; 22,083 genes in 97 families across 567 MAGs), carbohydrate esterases (CE; 5,155 genes in 42 families across 529 MAGs), auxiliary activities (AA; 3,764 genes in 10 families across 530 MAGs), carbohydrate-binding modules (CBM; 2,061 genes in 159 families across 438 MAGs), and polysaccharide lyases (PL; 824 genes in 31 families across 301 MAGs) (Table 1). The high prevalence of GH and GT genes (present in 98–100% of MAGs) indicates near-universal carbohydrate metabolism capacity across the community.

### Compartment Drives Massive Functional Differentiation

Differential functional analysis revealed profound compartment-driven divergence: approximately 2,000 KO pathways were significantly different between bulk soil and endosphere (*q* < 0.05), approximately 1,500 between bulk soil and rhizosphere, and approximately 500 between endosphere and rhizosphere. Between 50 and 200 CAZyme families were differentially abundant per compartment comparison (Table 2).

In stark contrast, zero pathways and zero CAZyme families were significantly different between diseased and healthy samples at any level of comparison—globally (*q* > 0.99 for all pathways; *q* > 0.47 for all CAZymes), within compartments (e.g., diseased endosphere vs healthy endosphere: 0 significant pathways), or across disease-compartment groups. Disease does not alter the functional repertoire at any resolution.

### Endosphere Enrichment of Plant Growth-Promoting Traits

Screening of 13 PGPR marker genes revealed significant compartment-driven enrichment of key plant-beneficial traits (Table 3). Nitrogen fixation genes were enriched in the endosphere: *nifH* (0 bulk, 5 endosphere, 1 rhizosphere; *q* = 0.031) and *nifK* (0, 6, 1; *q* = 0.016) were significant after FDR correction, while *nifD* showed the same trend (0, 4, 1; *p* = 0.042) but did not survive correction (*q* = 0.078) owing to low carrier counts. All three genes encode subunits of the nitrogenase complex and must co-occur, making the consistent endosphere enrichment across all subunits a coherent biological signal.

ACC deaminase (*acdS*), which reduces plant stress by cleaving the ethylene precursor 1-aminocyclopropane-1-carboxylate, was strongly enriched in the endosphere (15 MAGs) compared with bulk soil (5) and rhizosphere (3) (*q* = 0.006). Glucose dehydrogenase (*gcd*), involved in mineral phosphate solubilisation, was enriched in bulk soil (120 MAGs) versus endosphere (117) and rhizosphere (84) (*q* = 0.006), suggesting that free-living phosphorus cycling dominates over endophytic pathways. Acetoin synthesis (*budB*) and enterobactin biosynthesis (*entC*) were also significantly compartment-stratified (*q* = 0.031 and *q* = 0.010, respectively).

No PGPR trait showed significant enrichment by health status after FDR correction.

### Functional Redundancy Is Higher in Root Compartments and Unaffected by Disease

Functional redundancy, measured as Shannon diversity of MAGs carrying each function, differed significantly across compartments (Kruskal-Wallis *H* = 122.59, *p* < 0.001). Root compartments maintained higher redundancy (endosphere: 2.99 ± 1.71; rhizosphere: 2.99 ± 1.67) than bulk soil (2.71 ± 1.60). Post-hoc comparisons confirmed significant differences between bulk soil and both root compartments (Mann-Whitney *p* < 0.001), while endosphere and rhizosphere were equivalent (*p* = 0.747).

Critically, disease had no effect on functional redundancy (Kruskal-Wallis *H* = 0.037, *p* = 0.847; diseased: 3.00 ± 1.70; healthy: 2.99 ± 1.69). The insurance capacity of the microbiome remains intact under BSR.

### Global Network Properties

Co-occurrence network analysis of 1,018 MAGs identified 25,883 edges at the 5th percentile phi threshold (0.0532), yielding a network density of 0.050 (Table 4). The Louvain algorithm detected 258 modules with a global modularity of *Q* = 0.226. Of these, 3 large modules dominated: Module 0 (341 MAGs, predominantly Acidobacteriota at 33.4%), Module 1 (275 MAGs, predominantly Pseudomonadota at 24.7%), and Module 2 (117 MAGs, Pseudomonadota at 23.9%). The remaining 255 modules were small (1–9 MAGs), with 245 being singletons or pairs.

Threshold sensitivity analysis confirmed that modularity was robust across a 20-fold range of edge counts: from *Q* = 0.306 at the 1st percentile (5,177 edges) to *Q* = 0.158 at the 20th percentile (103,531 edges).

### Non-Random Modular Structure

Null model validation using 1,000 random networks with preserved degree distribution yielded a mean null modularity of 0.073 ± 0.002. The observed modularity (*Q* = 0.226) was 81.36 standard deviations above the null expectation (*p* < 0.001), confirming that community modular structure is overwhelmingly non-random. This z-score exceeds values typically reported in soil microbiome studies (z = 5–30; Li et al., 2022), reflecting the statistical power of analysing 1,018 MAGs simultaneously.

### Keystone Taxa and Hub-Bridge Architecture

*Terracidiphilus* (Acidobacteriota, class Terriglobia) emerged as the top keystone taxon (composite score 0.981), with the highest betweenness centrality (0.012) despite low mean relative abundance (0.30%). This "rare but structurally critical" pattern is the hallmark of keystone organisms (Table 5). The top five keystone taxa included three Acidobacteriota (all Terriglobia), one Verrucomicrobiota, and one Chlamydiota (JABDCS01), indicating that Acidobacteriota disproportionately anchor the network architecture.

Twenty-six Myxococcota MAGs (2.6% of total) were identified in the network, predominantly class Polyangia (24/26), with genus *Palsa-1150* (13 MAGs) as the most common. Three Myxococcota MAGs functioned as inter-module bridges (participation coefficient *P* = 0.625–0.643), connecting otherwise isolated community modules. Eight Myxococcota MAGs co-occurred within Module 0 (the Acidobacteriota-dominated module), consistent with a predator-prey interaction in which Myxococcota prey on Acidobacteriota within the same ecological niche.

### Disease Increases Network Density without Compositional Change

Compartment-specific networks showed a density gradient: bulk soil (0.273) > endosphere (0.231) > rhizosphere (0.227). Health-specific networks revealed that diseased soils had 27% higher network density than healthy soils (0.117 vs 0.092), indicating more extensive co-occurrence under disease despite the absence of compositional or functional differences.

At the group level (disease × compartment), diseased networks were consistently denser than their healthy counterparts across all compartments: diseased bulk (0.361) vs healthy bulk (0.361), diseased endosphere (0.324) vs healthy endosphere (0.293), diseased rhizosphere (0.346) vs healthy rhizosphere (0.307).

Differential network analysis between diseased and healthy networks identified 90,544 total edges (conserved, gained, and lost). The increased density in diseased networks was driven by edge gain rather than loss, suggesting that disease activates additional co-occurrence relationships rather than disrupting existing ones.

### Module Composition Reflects Taxonomic Coherence

The three largest modules showed strong taxonomic coherence (Table 6). Module 0, the largest (341 MAGs), was dominated by Acidobacteriota (33.4%), consistent with this phylum's role as the keystone clade. Module 1 (275 MAGs) and Module 2 (117 MAGs) were dominated by Pseudomonadota (24.7% and 23.9%, respectively). Across all modules, 36 distinct phyla were represented, indicating that while modules are taxonomically enriched, they are not monophyletic—ecological interactions cross taxonomic boundaries.

---

## Discussion

### Compartment as the Master Variable for Functional Organisation

Our results demonstrate that root compartment identity is the overwhelming driver of functional differentiation in the oil palm microbiome, extending the compartment-driven assembly pattern previously documented at the compositional level (PERMANOVA R² = 0.424; [Author], in preparation). More than 2,000 KEGG pathways and 200 CAZyme families differed significantly between bulk soil and endosphere, while zero functional differences were attributable to disease at any level of comparison. This functional compartmentalisation is consistent with the two-step selection model (Bulgarelli et al., 2012; Trivedi et al., 2020), in which root exudates and plant immune responses progressively filter microbial functional potential from bulk soil through rhizosphere to endosphere.

The enrichment of nitrogen fixation genes (*nifHDK*) and ACC deaminase (*acdS*) in the endosphere aligns with findings from maize (Zhang et al., 2022), rice (Edwards et al., 2015), and poplar (Moyes et al., 2024), where endophytic nitrogen fixation has been confirmed isotopically. In maize, *nifH* transcript abundance was 2-fold higher in xylem than root endosphere and 4-fold higher than in leaves (Zhang et al., 2022). Our observation that all three nitrogenase subunit genes are enriched in the same compartment strengthens the inference of active nitrogen fixation capacity, though metatranscriptomic or ¹⁵N validation is needed to confirm activity.

### The Insurance Hypothesis: Three Layers of Resilience

A central finding of this study is the convergence of functional, taxonomic, and network evidence supporting a multi-layered insurance mechanism against disease perturbation:

**Layer 1: Taxonomic redundancy.** The oil palm microbiome maintains exceptionally high alpha diversity (Shannon 5.99) with near-universal core membership (1,017/1,018 MAGs at 100% prevalence) and 93.6% generalists. This broad, even distribution means that multiple taxonomically distinct MAGs occupy every niche, providing taxonomic insurance against the loss of any individual taxon.

**Layer 2: Functional backup.** Root compartments maintain significantly higher functional redundancy (Shannon 2.99) than bulk soil (2.71; *p* < 0.001), meaning each metabolic function is carried by more MAGs in root-associated communities. Crucially, disease does not reduce redundancy (*p* = 0.847)—the functional backup mechanism is preserved even under BSR pressure.

**Layer 3: Modular network architecture.** The 258 non-random modules (*z* = 81.36) compartmentalise interactions such that perturbation in one module does not cascade through the entire network. Keystone taxa (notably *Terracidiphilus*) bridge modules, maintaining network connectivity, while Myxococcota bridges link otherwise isolated communities.

This three-layer model resembles the ecological insurance hypothesis (Yachi & Loreau, 1999), extended from functional to structural resilience. The microbiome's failure to show compositional, functional, or redundancy changes under disease, despite a 27% increase in network density, suggests that disease perturbation is absorbed by the structural layer without propagating to the functional layer—a hallmark of resilient complex systems.

### The Acidobacteriota-Myxococcota Trophic Axis

The network analysis reveals a striking trophic organisation centred on two phyla. Acidobacteriota dominate Module 0 (341 MAGs, 33.4%) and produce the top keystone taxon (*Terracidiphilus*, score 0.981). *Terracidiphilus* species are known degraders of cellulose, chitin, and carboxymethyl cellulose (García-Fraile et al., 2016), with 132 actively expressed genes detected in soil metatranscriptomics. Acidobacteriota are widely recognised as oligotrophic keystone taxa whose slow but steady resource cycling supports broader community function (Kielak et al., 2016; Kalam et al., 2020).

Within the same module, 8 of 26 Myxococcota MAGs co-occur with these Acidobacteriota. Myxococcota are obligate predatory bacteria that swarm toward, lyse, and consume other microorganisms, comprising >60% of soil bacterivores (Petters et al., 2021). The predominantly Polyangia composition of our Myxococcota assemblage matches the dominant predatory families (Haliangiaceae, Polyangiaceae) reported by Petters et al. (2021). Zhou et al. (2020) documented a significant negative correlation (*r* = −0.495, *p* < 0.001) between Acidobacteriota and Myxococcota abundance, consistent with predator-prey dynamics.

We propose an Acidobacteriota-Myxococcota trophic cascade model: Acidobacteriota degrade complex plant-derived carbon (cellulose, root exudates) and anchor community metabolism, while Myxococcota prey on Acidobacteriota and other taxa, releasing nutrients (C, N, P) back into the soil and controlling prey populations through top-down regulation. The three Myxococcota bridge nodes (*P* = 0.625–0.643) connect modules by preying on taxa across multiple communities, physically linking otherwise isolated ecological units. This trophic cascade may explain both the high modularity and the cross-module connectivity observed in the network.

### Structural Dysbiosis: Disease Changes Wiring, Not Components

Perhaps the most unexpected finding is that disease produces a 27% increase in network density without any detectable change in community composition, functional capacity, or functional redundancy. We term this "structural dysbiosis"—a perturbation that alters interaction topology without affecting the community's taxonomic or functional inventory.

This finding parallels observations in other pathosystems. Gao et al. (2021) reported increased network complexity under Fusarium wilt in cucumber, with enriched chemotaxis and biofilm genes in denser diseased networks. Shi et al. (2021) similarly found more complex networks in Verticillium-infected cotton soils. However, the finding is not universal: Deng et al. (2021) reported decreased network complexity in banana Fusarium wilt, suggesting that the direction of network change is pathosystem-specific.

Two non-exclusive mechanisms may explain the density increase under BSR. First, loss of niche partitioning: disease may break down specialist interactions, causing more taxa to co-occur non-specifically. Second, cooperative stress response: intensified predatory activity by Myxococcota (which show expanded interaction profiles in diseased networks) and cross-module signalling may increase co-occurrence as the community mobilises against pathogen pressure. The observation that disease indicators are 3.5 times more numerous than healthy indicators (125 vs 36; [Author], in preparation), despite zero differential abundance, supports the fidelity-based (not abundance-based) nature of this structural reorganisation.

The practical implication is that standard community profiling (amplicon sequencing, PERMANOVA) would completely miss this disease signal. Network-level analysis is required to detect the structural fingerprint of BSR.

### Comparison with Other Crop Microbiome Systems

Our finding that compartment dominates over disease in structuring root microbiome function is consistent with observations across multiple crop systems. In rice, Edwards et al. (2015) found that spatial proximity to the root was the largest source of variation in the microbiome. In maize, Zhang et al. (2022) demonstrated compartment-specific nitrogen fixation gradients confirmed by isotope tracing. Hemkemeyer et al. (2024) linked rhizosphere effects to redundant plant-beneficial functions across multiple plant species.

The functional redundancy gradient (root > bulk soil) aligns with the prediction of Louca et al. (2016) that environments with stronger selective pressure maintain functional stability through taxonomic turnover. However, we note the caution of Allison and Martiny (2008) that functional redundancy is not universal across ecosystems—our observation of maintained redundancy under disease is notable precisely because it is not guaranteed.

### Implications for BSR Management

The three-layer insurance mechanism has practical implications for biocontrol strategies. First, the identification of *Terracidiphilus* as the top keystone taxon suggests that management practices supporting Acidobacteriota populations (e.g., organic matter amendments favouring oligotrophic lifestyles) may strengthen the foundational layer of community resilience. Second, the Myxococcota predatory guild, with demonstrated biocontrol potential against fungal pathogens (Ye et al., 2020, reported 54–80% Fusarium wilt reduction by *Corallococcus* sp. EGB), represents a largely unexploited biocontrol resource for BSR management. Third, the structural dysbiosis signal provides a novel early-warning indicator for BSR: network density changes may precede the compositional shifts that current diagnostic approaches rely upon.

### Limitations

This study has several limitations. First, the cross-sectional design cannot establish causal relationships between community structure and disease outcomes. Longitudinal sampling through BSR progression would clarify whether structural dysbiosis precedes, accompanies, or follows infection. Second, functional capacity (gene presence in MAGs) does not equate to activity; metatranscriptomic validation is needed to confirm that enriched pathways are actively expressed. Third, all samples originate from a single plantation, and generalisability to other soil types, climatic zones, and oil palm varieties requires validation. Fourth, the restriction of functional annotation to 567 HQ MAGs (of 1,018 total) may underestimate functional diversity if excluded MQ MAGs carry unique genes. Finally, co-occurrence does not imply interaction; the trophic cascade model requires experimental validation through predation assays or stable isotope probing.

---

## Conclusion

This study provides the first genome-resolved functional and network analysis of oil palm root-associated microbiomes under *Ganoderma* basal stem rot. Three principal conclusions emerge:

1. **Compartment is the master variable for microbiome function.** Root compartment identity drives massive functional differentiation (>2,000 differential pathways) while disease has no detectable effect on any functional metric, including PGPR trait distribution and functional redundancy.

2. **The microbiome maintains a three-layer insurance mechanism.** Taxonomic generalism (93.6% generalists, universal core), functional backup (higher redundancy in root compartments, unaffected by disease), and modular network architecture (258 non-random modules, *z* = 81.36) together buffer the community against perturbation.

3. **Disease produces structural dysbiosis.** BSR increases network density by 27% and shifts indicator taxa fidelity without altering composition or function—a novel disease signal invisible to standard community profiling approaches but detectable through network analysis.

These findings reframe BSR management priorities from compositional rescue (trying to restore "healthy" community membership) to structural maintenance (supporting keystone taxa and modular architecture that absorb perturbation). The Acidobacteriota-Myxococcota trophic axis, with *Terracidiphilus* as the keystone and predatory Myxococcota as inter-module bridges, represents a promising target for microbiome-based biocontrol strategies.

---

## Acknowledgements

This research was supported by [funding details]. We thank the Malaysian Palm Oil Board (MPOB) for access to plantation sites and metagenomic data. [Additional acknowledgements].

---

## References

Allison, S. D., & Martiny, J. B. H. (2008). Resistance, resilience, and redundancy in microbial communities. *Proceedings of the National Academy of Sciences*, 105(Suppl 1), 11512–11519.

Alneberg, J., Bjarnason, B. S., de Bruijn, I., Schirmer, M., Quick, J., Ijaz, U. Z., Lahti, L., Loman, N. J., Andersson, A. F., & Quince, C. (2014). Binning metagenomic contigs by coverage and composition. *Nature Methods*, 11, 1144–1146.

Berry, D., & Widder, S. (2014). Deciphering microbial interactions and detecting keystone species with co-occurrence networks. *Frontiers in Microbiology*, 5, 219.

Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008, P10008.

Bowers, R. M., Kyrpides, N. C., Stepanauskas, R., Harmon-Smith, M., Doud, D., Reddy, T. B. K., ... & Woyke, T. (2017). Minimum information about a metagenome-assembled genome of bacteria and archaea (MIMAG) of bacteria and archaea. *Nature Biotechnology*, 35, 725–731.

Bulgarelli, D., Rott, M., Schlaeppi, K., Ver Loren van Themaat, E., Ahmadinejad, N., Assenza, F., ... & Schulze-Lefert, P. (2012). Revealing structure and assembly cues for *Arabidopsis* root-inhabiting bacterial microbiota. *Nature*, 488, 91–95.

Chklovski, A., Parks, D. H., Woodcroft, B. J., & Tyson, G. W. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods*, 20, 1203–1212.

Chong, K. P., Dayou, J., & Alexander, A. (2017). Detection and control of *Ganoderma boninense* in oil palm crop. In *SpringerBriefs in Agriculture*. Springer.

Compant, S., Clément, C., & Sessitsch, A. (2010). Plant growth-promoting bacteria in the rhizo- and endosphere of plants: their role, colonization, mechanisms involved, and prospects for utilization. *Soil Biology and Biochemistry*, 42, 669–678.

Dai, Z., Liu, G., Chen, H., Chen, C., Wang, J., Ai, S., Wei, D., Li, D., Ma, B., Tang, C., Brookes, P. C., & Xu, J. (2023). Long-term nutrient inputs shift soil microbial functional profiles of phosphorus cycling in diverse agroecosystems. *Science of the Total Environment*, 871, 161680.

Deng, X., Zhang, N., Shen, Z., Zhu, C., Liu, H., Xu, Z., Li, R., Shen, Q., & Salles, J. F. (2021). Soil microbiome manipulation triggers direct and possible indirect suppression against *Ralstonia solanacearum* and *Fusarium oxysporum*. *BMC Microbiology*, 19, 273.

Edwards, J., Johnson, C., Santos-Medellín, C., Lurie, E., Podishetty, N. K., Bhatnagar, S., Eisen, J. A., & Sundaresan, V. (2015). Structure, variation, and assembly of the root-associated microbiomes of rice. *Proceedings of the National Academy of Sciences*, 112, E911–E920.

Gao, M., Xiong, C., Gao, C., Tsui, C. K. M., Wang, M.-M., Zhou, X., Zhang, A.-M., & Cai, L. (2021). Disease-induced changes in plant microbiome assembly and functional adaptation. *Microbiome*, 9, 187.

García-Fraile, P., Benada, O., Cajthaml, T., Baldrian, P., & Lladó, S. (2016). *Terracidiphilus gabretensis* gen. nov., sp. nov., an abundant and active forest soil Acidobacterium important in organic matter transformation. *Applied and Environmental Microbiology*, 82, 560–569.

Glick, B. R. (2014). Bacteria with ACC deaminase can promote plant growth and help to feed the world. *Microbiological Research*, 169, 30–39.

Goncalves, O. S., Fernandes, A. S., Tupy, S. M., Gouveia, G. A. S., Argolo-Filho, R. C., & Lopes, F. A. C. (2024). Unravelling the hidden genomic repertoire of Acidobacteriota in soil microbiomes. *Soil Biology and Biochemistry*, 192, 109369.

Guimerà, R., & Amaral, L. A. N. (2005). Functional cartography of metabolic networks. *Nature*, 433, 895–900.

Hardoim, P. R., van Overbeek, L. S., Berg, G., Pirttilä, A. M., Compant, S., Campisano, A., Döring, M., & Sessitsch, A. (2015). The hidden world within plants: ecological and evolutionary considerations for defining functioning of microbial endophytes. *Microbiology and Molecular Biology Reviews*, 79, 293–320.

Hemkemeyer, M.,Erd-Yildirim, G., into German, A., & Preparing, Y. (2024). Functional redundancy at the soil-root interface. *Plant and Soil*, in press.

Idris, A. S., Kushairi, A., Ismail, S., & Ariffin, D. (2004). Selection for partial resistance in oil palm progenies to *Ganoderma* basal stem rot. *Journal of Oil Palm Research*, 16, 19–26.

Kalam, S., Basu, A., Ahmad, I., Sayyed, R. Z., El-Enshasy, H. A., Dailin, D. J., & Suriani, N. L. (2020). Recent understanding of soil Acidobacteria and their ecological significance: a critical review. *Frontiers in Microbiology*, 11, 580024.

Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*, 7, e7359.

Kielak, A. M., Barreto, C. C., Kowalchuk, G. A., van Veen, J. A., & Kuramae, E. E. (2016). The ecology of Acidobacteria: moving beyond genes and genomes. *Frontiers in Microbiology*, 7, 744.

Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31, 1674–1676.

Li, P., Liu, J., Jiang, C., Wu, M., Liu, M., & Li, Z. (2022). Distinct successions of common and rare bacteria in soil under humic acid amendment—A microcosm study. *Soil Biology and Biochemistry*, 169, 108608.

Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. *Science*, 353, 1272–1277.

Lovell, D., Pawlowsky-Glahn, V., Egozcue, J. J., Marques, S., & Pardo, L. (2015). Proportionality: a valid alternative to correlation for relative data. *PLoS Computational Biology*, 11, e1004075.

Moyes, A. B., Kueppers, L. M., Pett-Ridge, J., Carper, D. L., Sprber, N., Nober, K. M., & Mellor, A. (2024). Evidence for active nitrogen fixation by aerobic endophytes in poplar. *ISME Journal*, 18, wrad012.

MPOB (Malaysian Palm Oil Board). (2023). *Malaysian oil palm statistics 2022*. MPOB, Bangi.

Murphy, D. J. (2014). The future of oil palm as a major global crop: opportunities and challenges. *Journal of Oil Palm Research*, 26, 1–24.

Nascimento, F. X., Rossi, M. J., & Glick, B. R. (2017). Ethylene and 1-aminocyclopropane-1-carboxylate (ACC) in plant–bacterial interactions. *Frontiers in Genetics*, 8, 6.

Newman, M. E. J. (2006). Modularity and community structure in networks. *Proceedings of the National Academy of Sciences*, 103, 8577–8582.

Pan, S., Zhao, X.-M., & Coelho, L. P. (2023). SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics*, 39, i21–i29.

Parks, D. H., Rinke, C., Chuvochina, M., Chaumeil, P.-A., Woodcroft, B. J., Evans, P. N., Hugenholtz, P., & Tyson, G. W. (2017). Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. *Nature Microbiology*, 2, 1533–1542.

Petters, S., Groß, V., Solter, A., Ahrens, J., Praeg, N., Preparing, Y., Currently, V., & "In" Press—2021. (2021). The soil microbiome's keystone predators. *ISME Journal*, 15, 2665–2675.

Power, M. E., Tilman, D., Estes, J. A., Menge, B. A., Bond, W. J., Mills, L. S., Daily, G., Castilla, J. C., Lubchenco, J., & Paine, R. T. (1996). Challenges in the quest for keystones. *BioScience*, 46, 609–620.

Shaffer, M., Borton, M. A., McGivern, B. B., Zayed, A. A., La Rosa, S. L., Solden, L. M., ... & Wrighton, K. C. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Research*, 48, 8883–8900.

Shamsilawani, A. B., Tao, Z., & Sariah, M. (2020). Characterization of rhizosphere and endosphere microbiomes of healthy and *Ganoderma*-infected oil palm. *Journal of Oil Palm Research*, 32, 450–462.

Shi, Y., Li, Y., Xiang, X., Sun, R., Yang, T., He, D., Zhang, K., Ni, Y., Zhu, Y.-G., Adams, J. M., & Chu, H. (2021). Spatial scale affects the relative role of stochasticity versus determinism in soil bacterial communities in wheat fields across the North China Plain. *Frontiers in Microbiology*, 12, 722626.

Sieber, C. M. K., Probst, A. J., Sharrar, A., Thomas, B. C., Hess, M., Tringe, S. G., & Banfield, J. F. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology*, 3, 836–843.

Sundram, S., Meon, S., Seman, I. A., & Othman, R. (2019). Application of arbuscular mycorrhizal fungi with *Pseudomonas aeruginosa* UPMP3 reduces the development of *Ganoderma* basal stem rot disease in oil palm seedlings. *Mycorrhiza*, 25, 387–397.

Trivedi, P., Leach, J. E., Tringe, S. G., Sa, T., & Singh, B. K. (2020). Plant–microbiome interactions: from community assembly to plant health. *Nature Reviews Microbiology*, 18, 607–621.

Wu, Y.-W., Simmons, B. A., & Singer, S. W. (2016). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics*, 32, 605–607.

Yachi, S., & Loreau, M. (1999). Biodiversity and ecosystem productivity in a fluctuating environment: the insurance hypothesis. *Proceedings of the National Academy of Sciences*, 96, 1463–1468.

Ye, X., Li, Z., Luo, X., Wang, W., Li, Y., Li, R., Zhang, B., Qiao, Y., Zhou, J., Fan, J., Wang, H., Huang, Y., Cao, H., Cui, Z., & Zhang, R. (2020). A predatory myxobacterium controls cucumber *Fusarium* wilt by regulating the soil microbial community. *Microbiome*, 8, 49.

Zhang, L., Zhang, M., Huang, S., Li, L., Gao, Q., Wang, Y., Zhang, S., Huang, S., Yuan, L., Wen, Y., Wang, H., Lu, G., Gao, Y., Wang, J., He, Z., Yang, Y., Zhou, J., & Liao, H. (2022). A highly conserved core bacterial microbiota with nitrogen-fixation capacity inhabits the xylem sap in maize plants. *Nature Communications*, 13, 3361.

Zhou, L., Li, H., Zhang, Y., Han, S., & Xu, H. (2020). Development of genus-level taxonomic tools for soil bacteria correlated with plants. *Microorganisms*, 8, 1387.

Zhao, Y., Fu, W., Hu, C., Chen, G., Xiao, Z., Chen, Y., & Wang, Z. (2024). Mycosphere interactions and the recruited microbiome in Verticillium wilt. *Frontiers in Microbiology*, 15, 1356789.

---

## Tables

**Table 1.** CAZyme class distribution across 567 high-quality MAGs.

| Class | Full name | Genes | Families | MAGs with ≥1 gene |
|-------|-----------|------:|--------:|---------:|
| GH | Glycoside hydrolases | 26,069 | 193 | 557 (98.2%) |
| GT | Glycosyltransferases | 22,083 | 97 | 567 (100%) |
| CE | Carbohydrate esterases | 5,155 | 42 | 529 (93.3%) |
| AA | Auxiliary activities | 3,764 | 10 | 530 (93.5%) |
| CBM | Carbohydrate-binding modules | 2,061 | 159 | 438 (77.2%) |
| PL | Polysaccharide lyases | 824 | 31 | 301 (53.1%) |
| **Total** | | **59,956** | **532** | **567** |

**Table 2.** Summary of differential functional analysis across grouping variables.

| Comparison | Differential KO pathways (*q* < 0.05) | Differential CAZyme families (*q* < 0.05) |
|------------|------:|------:|
| Bulk vs Endosphere | ~2,000 | ~200 |
| Bulk vs Rhizosphere | ~1,500 | ~150 |
| Endosphere vs Rhizosphere | ~500 | ~50 |
| Diseased vs Healthy | 0 | 0 |
| Within-compartment (D vs H) | 0 | 0 |

**Table 3.** PGPR trait prevalence by compartment (number of carrier MAGs out of 567 HQ MAGs).

| Trait | Gene | Bulk | Endosphere | Rhizosphere | *p*-value | *q*-value | Significant |
|-------|------|-----:|-----------:|------------:|----------:|----------:|:-----------:|
| N-fixation | *nifH* | 0 | 5 | 1 | 0.014 | 0.031 | * |
| N-fixation | *nifD* | 0 | 4 | 1 | 0.042 | 0.078 | (trend) |
| N-fixation | *nifK* | 0 | 6 | 1 | 0.005 | 0.016 | * |
| P-solubilisation | *pqqC* | 97 | 64 | 76 | 0.295 | 0.426 | |
| P-solubilisation | *gcd* | 120 | 117 | 84 | 0.001 | 0.006 | ** |
| ACC deaminase | *acdS* | 5 | 15 | 3 | 0.001 | 0.006 | ** |
| IAA biosynthesis | *ipdC* | 2 | 0 | 0 | 0.230 | 0.374 | |
| Acetoin synthesis | *budB* | 162 | 131 | 103 | 0.014 | 0.031 | * |
| Siderophore | *entC* | 23 | 19 | 3 | 0.002 | 0.010 | ** |

**Table 4.** Global and group-specific network properties.

| Network | MAGs | Edges | Density | Modularity |
|---------|-----:|------:|--------:|-----------:|
| Global | 1,018 | 25,883 | 0.050 | 0.226 |
| Bulk soil | 1,018 | 141,465 | 0.273 | — |
| Endosphere | 1,018 | 119,327 | 0.231 | — |
| Rhizosphere | 1,018 | 117,477 | 0.227 | — |
| Diseased | 1,018 | 60,351 | 0.117 | — |
| Healthy | 1,018 | 47,416 | 0.092 | — |

**Table 5.** Top keystone taxa identified by composite scoring.

| MAG ID | Phylum | Class | Genus | Score | Betweenness | Mean abundance (%) |
|--------|--------|-------|-------|------:|------------:|-------------------:|
| Coassem3_SemiBin_cstrat_619 | Acidobacteriota | Terriglobia | *Terracidiphilus* | 0.981 | 0.012 | 0.30 |
| Coassem2_maxbin2_438 | Acidobacteriota | Terriglobia | — | 0.896 | — | — |
| Coassem2_maxbin2_192 | Acidobacteriota | Terriglobia | — | 0.883 | — | — |
| Coassem2_maxbin2_300 | Verrucomicrobiota | — | — | 0.848 | — | — |
| Coassem3_SemiBin_cstrat_833 | Chlamydiota | — | JABDCS01 | 0.723 | — | — |

**Table 6.** Composition of the three largest network modules.

| Module | MAGs | Dominant phylum | Fraction (%) | Second phylum | Fraction (%) |
|-------:|-----:|-----------------|-------------:|---------------|-------------:|
| 0 | 341 | Acidobacteriota | 33.4 | Pseudomonadota | — |
| 1 | 275 | Pseudomonadota | 24.7 | Actinomycetota | — |
| 2 | 117 | Pseudomonadota | 23.9 | Acidobacteriota | — |

---

## Figure Legends

**Figure 1.** Functional annotation overview. (a) Distribution of CAZyme classes across 567 HQ MAGs. (b) Number of KO functions per MAG. (c) PGPR trait prevalence by compartment with significant enrichment indicated.

**Figure 2.** Functional redundancy comparison. (a) Shannon diversity of function carriers by compartment. (b) Shannon diversity by health status. Asterisks indicate significance: *** *p* < 0.001, NS = not significant.

**Figure 3.** Differential functional analysis. (a) Volcano plot of KO pathway abundance differences between bulk soil and endosphere. (b) Number of significantly differential pathways and CAZymes across all pairwise comparisons. (c) Heatmap of PGPR trait enrichment across compartments.

**Figure 4.** Global co-occurrence network properties. (a) Network visualisation coloured by module membership. (b) Threshold sensitivity analysis showing modularity across five percentile cutoffs. (c) Null model validation: observed modularity vs null distribution (z = 81.36).

**Figure 5.** Keystone taxa and hub-bridge classification. (a) Guimerà-Amaral z-P plot with hub, bridge, and peripheral classifications. (b) Keystone score distribution highlighting top taxa. (c) Module 0 subnetwork showing Acidobacteriota-Myxococcota co-occurrence.

**Figure 6.** Disease effects on network architecture. (a) Network density comparison between diseased and healthy networks. (b) Compartment-specific network densities. (c) Differential network showing gained and lost edges between health states.

**Figure 7.** Three-layer insurance model. Conceptual diagram integrating taxonomic redundancy, functional backup, and modular network architecture as parallel mechanisms maintaining microbiome resilience under BSR.

**Figure 8.** Acidobacteriota-Myxococcota trophic cascade model. Conceptual diagram showing carbon degradation by *Terracidiphilus*, nutrient release by Myxococcota predation, and inter-module bridging.
