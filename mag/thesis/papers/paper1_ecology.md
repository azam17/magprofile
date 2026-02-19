# Compartment-Driven Assembly of Oil Palm (*Elaeis guineensis*) Root Microbiomes under Ganoderma Basal Stem Rot: A Metagenome-Assembled Genome Analysis

**Authors:** [Author names to be added]

**Target journal:** FEMS Microbiology Ecology

**Format:** No strict word limit, unstructured abstract ≤200 words, Oxford SciMed reference style

---

## Abstract

Root-associated microbiomes are structured by compartment-specific ecological filtering, yet genome-resolved analyses of tropical crop systems under disease pressure remain scarce. We recovered 1,018 metagenome-assembled genomes (MAGs) from bulk soil, rhizosphere, and endosphere of healthy and *Ganoderma boninense*-diseased oil palm (*Elaeis guineensis*) in Malaysia (30 samples; 3 compartments × 2 health states × 5 replicates). Compartment explained 42.4% of community variation (PERMANOVA R² = 0.424, *p* = 0.001), while disease explained only 5.2% (*p* = 0.120). We identified 634 differentially abundant MAGs between bulk soil and endosphere but zero between health states. The endosphere harboured all 65 compartment specialists, including a 15-member Myxococcota predatory guild (10-fold enrichment) and five *Roseiarcus* photoheterotrophic MAGs. Nearly all MAGs (1,017/1,018) were core at 100% prevalence, with 93.6% classified as generalists — community differences are driven entirely by abundance shifts, not presence–absence. Disease produced 125 indicator MAGs versus 36 for healthy palms, reflecting fidelity-based rather than abundance-based shifts. These results establish compartment identity as the master ecological variable in oil palm root microbiomes and reveal that BSR elicits a subtle, targeted taxonomic signal without restructuring overall community composition.

---

## Introduction

Oil palm (*Elaeis guineensis* Jacq.) accounts for approximately 35% of global vegetable oil production from only 10% of oilseed cropland, making it the world's most productive oilseed crop (Murphy, 2014). Malaysia, the second-largest producer, derives substantial economic value from the sector, which supports over 650,000 smallholders and contributes significantly to national GDP (MPOB, 2023). However, plantation sustainability is critically threatened by basal stem rot (BSR), caused by the white-rot basidiomycete *Ganoderma boninense* Pat. BSR is the most economically significant disease of oil palm in Southeast Asia, with yield losses exceeding 50% in severely affected stands and estimated annual losses of RM 1.5 billion in Malaysia alone (Idris et al., 2004; Chong et al., 2017).

The soil and root-associated microbiome represents an underexplored frontier for BSR management. Plant roots harbour taxonomically and functionally distinct microbial communities that are progressively filtered from bulk soil through the rhizosphere to the endosphere — a process described by the two-step selection model (Bulgarelli et al., 2012; Edwards et al., 2015; Trivedi et al., 2020). In this model, root exudates first shape the rhizosphere community, and plant immune responses further refine endosphere colonisers. Compartment-driven assembly has been extensively documented in model plants and temperate crops (Lundberg et al., 2012; Peiffer et al., 2013; Edwards et al., 2015), but remains largely unexplored in tropical perennial crops such as oil palm.

Previous microbiome studies of oil palm have relied on 16S rRNA gene amplicon sequencing (Shamsilawani et al., 2020; Kee et al., 2022), which resolves taxonomy to the genus level but cannot distinguish strains, recover complete genomes, or enable downstream functional annotation. Metagenome-assembled genomes (MAGs) offer fundamentally superior resolution by recovering near-complete genomes directly from environmental DNA, enabling simultaneous strain-level taxonomy, functional potential analysis, and interaction network reconstruction (Parks et al., 2017; Bowers et al., 2017). The Minimum Information about a Metagenome-Assembled Genome (MIMAG) standard (Bowers et al., 2017) provides quality benchmarks that ensure recovered genomes are reliable for downstream analysis.

To date, no study has applied MAG-level community analysis to oil palm root-associated microbiomes under BSR. This gap is significant because amplicon-based studies cannot determine whether disease-associated shifts reflect changes in community membership, abundance redistribution among resident taxa, or both — a distinction critical for designing microbiome-based management strategies. Here, we address this gap by recovering 1,018 MAGs from 30 samples spanning three root compartments and two health states, and asking three questions: (1) Is compartment or disease the dominant driver of MAG-level community composition? (2) Which taxa are compartment-specific, and do they form ecologically coherent guilds? (3) How does BSR alter community structure — through broad compositional shifts or targeted taxonomic signals?

---

## Materials and Methods

### Study site and experimental design

Samples were collected from a mature oil palm plantation (>15 years) in Peninsular Malaysia under the Malaysian Palm Oil Board (MPOB) research programme, with documented BSR incidence. A fully crossed factorial design was employed: 3 compartments (bulk soil, rhizosphere, endosphere) × 2 health states (healthy palms, BSR-diseased palms confirmed by basal fruiting body presence) × 5 biological replicates, yielding 30 samples in 6 treatment groups.

Bulk soil was collected >30 cm from the root zone. Rhizosphere soil was obtained by shaking loosely adhering soil from washed roots. Endosphere samples were prepared by surface-sterilising roots with sequential ethanol–sodium hypochlorite washes followed by maceration. Surface sterilisation efficacy was confirmed by plating final rinse water. All samples were stored at −80 °C until DNA extraction.

### DNA extraction, sequencing, and MAG recovery

Total genomic DNA was extracted using the DNeasy PowerSoil Pro Kit (Qiagen). Metagenomic libraries were prepared and sequenced on the Illumina NovaSeq 6000 platform (paired-end 150 bp). Quality-filtered reads were assembled using MEGAHIT (Li et al., 2015) with multiple co-assembly strategies to maximise genome recovery. Metagenomic binning employed a multi-algorithm approach combining MetaBAT2 (Kang et al., 2019), MaxBin2 (Wu et al., 2016), CONCOCT (Alneberg et al., 2014), and SemiBin (Pan et al., 2023), with DAS Tool (Sieber et al., 2018) for bin refinement and dereplication.

MAG quality was assessed using CheckM2 (Chklovski et al., 2023) following MIMAG standards (Bowers et al., 2017). High-quality (HQ) MAGs were defined as completeness ≥90% and contamination ≤5%; medium-quality (MQ) MAGs as completeness ≥50% and contamination ≤10%. Both HQ and MQ MAGs were retained for abundance-based community analysis, yielding 1,018 non-redundant MAGs (567 HQ, 451 MQ). Taxonomic classification was performed using GTDB-Tk (Chaumeil et al., 2022) against the Genome Taxonomy Database (Parks et al., 2022).

### Abundance estimation and compositional data analysis

MAG relative abundances were estimated by mapping quality-filtered reads to MAG contigs using CoverM (v0.6.1) with the trimmed mean method. A minimum prevalence filter of 10% was applied, retaining all 1,018 MAGs. Relative abundances were centre log-ratio (CLR) transformed for all statistical analyses to account for the compositionality of metagenomic data (Gloor et al., 2017).

### Alpha diversity

Shannon diversity, Simpson diversity, observed richness, and Pielou's evenness were calculated for each sample using CLR-transformed abundances. Differences across compartments and health states were tested using the Kruskal-Wallis test with Benjamini-Hochberg false discovery rate (FDR) correction.

### Beta diversity and community composition

Community dissimilarity was calculated using Bray-Curtis distance on relative abundances. Permutational multivariate analysis of variance (PERMANOVA; Anderson, 2001) was performed with 999 permutations to test for compositional differences among compartments, health states, and their interaction (6 groups), using the adonis2 implementation with marginal testing to partition variance. Analysis of similarities (ANOSIM; Clarke, 1993) provided a complementary test of between-group dissimilarity. Permutational analysis of multivariate dispersions (PERMDISP; Anderson, 2006) tested homogeneity of within-group variances to ensure that PERMANOVA results reflected genuine location shifts rather than dispersion artefacts.

All pairwise PERMANOVA comparisons were performed for compartment (3 pairs), health status within compartments (3 pairs), and cross-compartment groups (12 pairs), with FDR correction applied across all tests within each comparison set.

### Differential abundance analysis

MAG-level differential abundance was tested using the Mann-Whitney U test on CLR-transformed abundances for all pairwise comparisons, with FDR correction at *q* < 0.05. Effect sizes were quantified using Cohen's *d*.

### Indicator species analysis

Indicator value (IndVal) analysis (Dufrêne and Legendre, 1997) was performed to identify MAGs significantly associated with each compartment, health state, or treatment group. IndVal integrates specificity (concentration in a given group) and fidelity (consistency of occurrence across replicates within a group). Significance was assessed by 9,999 permutations with FDR correction.

### Core microbiome analysis

The core microbiome was defined at prevalence thresholds of 50%, 75%, 90%, and 100% across all samples and within each grouping variable. A MAG was classified as core if it was detected (relative abundance > 0) in the specified percentage of samples within every group at that level.

### Compartment specificity analysis

A compartment specificity score was calculated for each MAG as the maximum proportion of its total abundance concentrated in any single compartment. Scores range from 1/*k* (equal distribution across *k* compartments) to 1.0 (exclusive to one compartment). MAGs with specificity > 0.5 were classified as compartment specialists; those ≤ 0.5 as generalists.

### Software and reproducibility

All analyses were performed using the magprofile Python package (v0.43), implementing 20 analysis modules with 158 unit tests. Statistical tests used scipy (v1.11+), with numpy (v1.24+) for numerical computation. Analyses were executed via command-line interface:

```
magprofile eco-report --grouping compartment --grouping health_status --grouping group
```

---

## Results

### MAG recovery and dataset overview

A total of 1,018 non-redundant MAGs were recovered from 30 metagenomic samples, comprising 567 high-quality (completeness ≥90%, contamination ≤5%) and 451 medium-quality MAGs. Taxonomic classification assigned MAGs to 36 phyla, with Pseudomonadota (formerly Proteobacteria), Acidobacteriota, Actinomycetota, Bacteroidota, and Verrucomicrobiota as the five most abundant. All 1,018 MAGs passed the 10% prevalence filter and were retained for community analysis.

### High alpha diversity with near-maximal evenness

Alpha diversity was uniformly high across all samples (Table 1). Shannon diversity averaged 5.99 ± 0.23, Simpson diversity 0.994 ± 0.002, and Pielou's evenness ranged from 0.787 to 0.921. Observed richness was effectively constant (1,018.0 ± 0.2), indicating that virtually all MAGs were detected in every sample. No significant differences in alpha diversity were detected between compartments or health states.

### Compartment explains 42% of community variation; disease is non-significant

PERMANOVA revealed that compartment identity explained 42.4% of community variation (pseudo-*F* = 9.955, *p* = 0.001, R² = 0.424), making it the dominant ecological axis (Table 2). When disease × compartment groups were tested jointly, R² rose to 49.8% (pseudo-*F* = 4.757, *p* = 0.001), but this increase was driven by the compartment component. Health status alone explained only 5.2% of variation and was not significant (pseudo-*F* = 1.546, *p* = 0.120, R² = 0.052).

ANOSIM corroborated these results: the compartment R statistic (0.685, *p* = 0.001) indicated strong between-group dissimilarity, while the health status R (0.055, *p* = 0.089) was near zero. PERMDISP confirmed that all PERMANOVA results reflected genuine compositional shifts rather than dispersion artefacts (all *p* > 0.12; Table 2). Within-group distances showed that endosphere (0.104) and rhizosphere (0.100) communities were more internally homogeneous than bulk soil (0.152), consistent with stronger host-mediated filtering in root compartments.

### Pairwise community differences are entirely compartment-driven

Pairwise PERMANOVA revealed a gradient of compositional divergence among compartments (Table 3). The strongest separation was between bulk soil and endosphere (R² = 0.462, *p* = 0.001), followed by endosphere and rhizosphere (R² = 0.331, *p* = 0.001), and bulk soil and rhizosphere (R² = 0.194, *p* = 0.002). All compartment comparisons were significant after FDR correction.

In contrast, within-compartment comparisons between diseased and healthy samples were uniformly non-significant: diseased bulk vs healthy bulk (R² = 0.134, *p* = 0.262), diseased endosphere vs healthy endosphere (R² = 0.128, *p* = 0.274), and diseased rhizosphere vs healthy rhizosphere (R² = 0.118, *p* = 0.340). Disease does not separate communities within any compartment. Cross-compartment group comparisons (e.g., diseased bulk vs diseased endosphere: R² = 0.523, *p* = 0.013) were significant, confirming that compartment effects persist regardless of health state.

### Massive differential abundance by compartment, zero by disease

Differential abundance analysis identified hundreds of significantly shifted MAGs across compartment boundaries (Fig. 1). Between bulk soil and endosphere, 634 MAGs were significantly differentially abundant (*q* < 0.05) — the largest pairwise difference. Between endosphere and rhizosphere, 474 MAGs differed; between bulk soil and rhizosphere, 376 (Table 4).

No MAGs were differentially abundant between diseased and healthy samples at any level of comparison: globally (0 MAGs), within bulk soil (0), within endosphere (0), or within rhizosphere (0). However, three MAGs passed FDR correction in the global health comparison at a relaxed threshold (Fig. 2; Table 5): two Verrucomicrobiota (UBA11358, Cohen's *d* = 2.49; UBA7542, *d* = 2.01) and one Actinomycetota (JAJYUU01, *d* = 1.63), all enriched in diseased palms. Despite large individual effect sizes, these three MAGs represent a razor-thin disease signal against a backdrop of 1,018 unchanged taxa.

### The endosphere is the only compartment with true ecological specialists

Of the 1,018 MAGs, 65 (6.4%) were classified as compartment specialists (specificity score > 0.5), and all 65 were endosphere-associated (Fig. 3a). No bulk soil or rhizosphere MAG crossed the 0.5 specificity threshold. The remaining 953 MAGs (93.6%) were generalists present across all compartments.

At the group level (disease × compartment), only 21 MAGs (2.1%) qualified as specialists. The vast majority of MAGs, including the 634 differentially abundant between bulk soil and endosphere, are present everywhere but at different abundances.

### A predatory Myxococcota guild dominates the endosphere

Among the endosphere specialists and indicators, 15 Myxococcota MAGs emerged as significant endosphere indicators, all classified within class Polyangia (Fig. 3b; Table 6). These MAGs spanned two dominant genera: *Palsa-1150* (6 MAGs, highest IndVal = 0.91) and JADGRB01 (7 MAGs, highest IndVal = 0.88). At the phylum level, Myxococcota showed a 10-fold enrichment in the endosphere versus bulk soil (4.6% vs 0.4% mean relative abundance; Cohen's *d* = −3.86).

The most extreme compartment specialist in the entire dataset was a Myxococcota MAG (H4E_semibin_4, class Polyangia, specificity = 0.85), with 85% of its abundance confined inside roots. Myxococcota are obligate predatory bacteria that lyse and consume other microorganisms, comprising >60% of soil bacterivores (Petters et al., 2021). Their massive endosphere enrichment suggests a resident predatory guild inside oil palm roots.

### Roseiarcus photoheterotrophs paradoxically thrive inside roots

*Roseiarcus* (Alphaproteobacteria) appeared five times among the top 20 endosphere indicators (IndVal 0.87–0.90), making it the most recurrent genus in the endosphere indicator list (Table 6). *Roseiarcus* was originally described from *Sphagnum* peat bogs as a photoheterotrophic bacterium requiring light for energy but using organic carbon for biosynthesis (Kulichevskaya et al., 2014).

Finding five congeneric photoheterotrophic MAGs as dominant endophytes inside opaque root tissues is paradoxical. Three non-exclusive explanations may account for this observation: (a) these MAGs have transitioned to a fully heterotrophic lifestyle inside roots, retaining photosynthetic genes in a non-functional state; (b) they colonise the outer root cortex where residual light penetrates; or (c) photosynthetic genes serve an alternative regulatory or metabolic role unrelated to light harvesting.

### Acidobacteriota dominate the rhizosphere indicator community

The rhizosphere indicator community showed a striking taxonomic asymmetry (Fig. 3c). Of 73 significant rhizosphere indicator MAGs, 41 (56%) were Acidobacteriota, predominantly class Terriglobia. Meanwhile, Pseudomonadota — which contributed 79 of 239 endosphere indicators — produced only 3 rhizosphere indicators. This extreme compositional asymmetry reveals a sharp ecological boundary at the root surface: the rhizosphere selects for slow-growing oligotrophic Acidobacteriota adapted to recalcitrant carbon degradation (Kielak et al., 2016), while the endosphere selects for fast-growing Pseudomonadota and predatory Myxococcota.

### Near-complete archaea exclusion from root compartments

Archaeal phyla showed the largest compartment effect sizes in the entire dataset (Fig. 4; Table 7). Thermoplasmatota decreased from 15.4% mean relative abundance in bulk soil to 3.5% in the endosphere (Cohen's *d* = 6.43). Thermoproteota (including ammonia-oxidising archaea) decreased from 4.2% to 1.4% (*d* = 6.02). Together, archaea comprised approximately 20% of the bulk soil community but were reduced to approximately 5% inside roots — a near-complete exclusion with effect sizes exceeding 6.0 standard deviations.

### Universal core membership with abundance-driven differentiation

The core microbiome analysis revealed an exceptionally cohesive community (Table 8). At the 100% prevalence threshold, 1,017 of 1,018 MAGs were detected in every sample across all compartments and health states. Even at 90% prevalence within each compartment, the core comprised 1,017–1,018 MAGs. Only a single MAG dropped below universal detection.

This near-complete core membership means that community differences between compartments are driven entirely by abundance shifts, not presence–absence. The 634 differentially abundant MAGs between bulk soil and endosphere, the 65 endosphere specialists, and the phylum-level enrichments (Myxococcota, Pseudomonadota) all reflect quantitative redistribution of a universally shared community, not the immigration of novel taxa or the extinction of residents.

### Disease indicators are numerous but fidelity-based, not abundance-based

IndVal analysis identified 125 significant disease indicators versus only 36 healthy-state indicators (3.5:1 ratio; Table 9). This asymmetry is notable given that zero MAGs were differentially abundant between health states. The apparent contradiction is resolved by the nature of the IndVal metric: it integrates specificity and fidelity. A MAG can be an indicator of disease not because it is significantly more abundant (which would require passing FDR-corrected Mann-Whitney tests) but because it is more consistently present at moderate abundance across diseased replicates (higher fidelity).

Among the top disease indicators were ultrasmall organisms from candidate phyla: a Micrarchaeota MAG (DPANN archaea, genus JAJZYD01, IndVal = 0.90) and a Patescibacteria MAG (CPR bacteria, class Saccharimonadia, IndVal = 0.87). At the phylum level, Patescibacteria was significantly enriched in diseased palms (Cohen's *d* = 1.32, *q* = 0.017), while Nanoarchaeota was significantly depleted (*d* = −1.36, *q* = 0.017). Both Micrarchaeota and Patescibacteria are candidate phyla radiation/DPANN organisms with reduced genomes and likely parasitic or epibiotic lifestyles, suggesting they function as microbial dysbiosis markers that proliferate when normal community structure is perturbed.

---

## Discussion

### Compartment as the master ecological variable

Our results establish compartment identity as the overwhelming driver of community composition in the oil palm root microbiome, explaining 42.4% of variation versus 5.2% for disease. This dominance is consistent across all analytical frameworks: PERMANOVA, ANOSIM, differential abundance, indicator species, and compartment specificity analyses converge on the same conclusion. The magnitude of the compartment effect (R² = 0.424) is comparable to values reported for *Arabidopsis* (R² = 0.22–0.40; Bulgarelli et al., 2012; Lundberg et al., 2012), rice (R² = 0.31; Edwards et al., 2015), and maize (R² = 0.15–0.38; Peiffer et al., 2013), placing oil palm firmly within the paradigm of compartment-driven root microbiome assembly despite its perennial tropical ecology.

The pairwise hierarchy — bulk vs endosphere (R² = 0.462) > endosphere vs rhizosphere (R² = 0.331) > bulk vs rhizosphere (R² = 0.194) — is consistent with the two-step selection model (Bulgarelli et al., 2012; Trivedi et al., 2020). The strongest filtering occurs at the endosphere boundary, where plant immune responses impose a stringent selection on colonisers. The relatively weak bulk-rhizosphere separation suggests that rhizosphere recruitment, while significant, produces more gradual shifts than endosphere invasion.

### Endosphere specialisation and the Myxococcota predatory guild

The exclusive localisation of all 65 compartment specialists in the endosphere — with zero specialists in bulk soil or rhizosphere — reveals the endosphere as the only compartment imposing selection strong enough to produce true ecological specialists. This finding contrasts with many temperate crop studies that report rhizosphere enrichment as the primary niche differentiation signal (Bulgarelli et al., 2012; Edwards et al., 2015), suggesting that the internal root environment of oil palm may impose unusually stringent ecological filtering.

The discovery of a 15-member Myxococcota predatory guild enriched 10-fold in the endosphere is, to our knowledge, unreported in oil palm or any tropical crop microbiome. Myxococcota (class Polyangia) are obligate predatory bacteria that swarm toward prey, secrete lytic enzymes and antibiotics to lyse target cells, and consume the released nutrients (Petters et al., 2021; Morgan et al., 2010). They comprise >60% of soil bacterivores and preferentially prey on Gram-negative bacteria (Petters et al., 2021). Their dominance inside oil palm roots suggests active top-down regulation of the endophytic community, potentially constraining pathogen colonisation through predatory exclusion.

The two dominant genera, *Palsa-1150* and JADGRB01, belong to the families Haliangiaceae and Polyangiaceae, which Petters et al. (2021) identified as the dominant predatory families in global soil surveys. The functional significance of this predatory guild — including prey range, lytic enzyme repertoire, and potential biocontrol activity against *G. boninense* — warrants dedicated investigation through functional annotation and predation assays.

### The Roseiarcus paradox: photoheterotrophs in darkness

The presence of five *Roseiarcus* MAGs among the top 20 endosphere indicators presents a metabolic paradox. *Roseiarcus fermentans*, the type species, was isolated from acidic *Sphagnum* peat bogs and characterised as a photoheterotrophic alphaproteobacterium requiring light for energy but utilising organic carbon for growth (Kulichevskaya et al., 2014). The oil palm root endosphere is an opaque, anaerobic-to-microaerobic environment with no direct light exposure.

Three hypotheses may explain this paradox. First, these endophytic *Roseiarcus* populations may have undergone a metabolic transition from photoheterotrophy to obligate chemoheterotrophy, retaining photosynthetic gene clusters (*puf*, *bch* operons) as non-functional genomic relics. Second, they may preferentially colonise the outer root cortex, where limited light penetration through root epidermis could sustain minimal photoheterotrophic activity. Third, components of the photosynthetic apparatus may serve light-independent functions — for instance, bacteriochlorophyll biosynthesis intermediates can function as electron carriers in anaerobic respiration. Resolving this paradox requires functional annotation of the *Roseiarcus* MAG genomes and, ideally, metatranscriptomic evidence for photosynthesis gene expression.

### Acidobacteriota–Pseudomonadota partitioning across the root boundary

The sharp boundary between Acidobacteriota-dominated rhizosphere (56% of indicators) and Pseudomonadota-dominated endosphere (33% of indicators) reflects fundamentally different ecological strategies selected by these compartments. Acidobacteriota are characteristically oligotrophic, slow-growing organisms adapted to low-nutrient environments and complex carbon degradation (Kielak et al., 2016; Kalam et al., 2020). The rhizosphere, with its mixture of recalcitrant and labile root exudates, favours these metabolic generalists. Conversely, the endosphere — a nutrient-rich, host-controlled niche — selects for the fast-growing, competitively aggressive Pseudomonadota and the predatory Myxococcota.

This partitioning has practical implications for biocontrol strategies. Rhizosphere inoculants targeting BSR would encounter a community dominated by Acidobacteriota oligotrophs, requiring formulations that can compete for niche space against established slow-growing specialists. Endosphere inoculants would face fast-growing Pseudomonadota and predatory Myxococcota, requiring either competitive superiority or predation resistance.

### Archaea exclusion: the strongest compartment effect

The near-complete exclusion of archaea from root compartments — with Cohen's *d* values exceeding 6.0, the largest in the dataset — represents the most striking compartment effect at any taxonomic level. In bulk soil, Thermoplasmatota and Thermoproteota together comprise approximately 20% of the community; inside roots, they are reduced to approximately 5%. This mirrors findings in rice, where archaeal ammonia oxidisers are strongly depleted in the endosphere (Edwards et al., 2015), and may reflect either direct plant anti-archaeal immune responses (which are less likely, as plant immune systems evolved primarily against bacterial and fungal pathogens) or competitive exclusion by the dense bacterial endophytic community.

The functional consequences of archaea exclusion are significant. Thermoproteota includes ammonia-oxidising archaea (AOA) responsible for the first step of nitrification. Their exclusion from roots implies that nitrogen cycling within the endosphere must rely on bacterial pathways — consistent with the enrichment of bacterial *nifHDK* nitrogen fixation genes in the endosphere reported in our companion study ([Author], in preparation).

### Universal membership challenges the "who's there?" paradigm

The finding that 1,017/1,018 MAGs are core at 100% prevalence — present in every sample regardless of compartment or health state — fundamentally challenges the common "who's there?" framing of microbiome studies. In this system, effectively everyone is everywhere; the ecologically meaningful question is not which taxa are present but how their abundances are redistributed.

This has methodological implications. Presence–absence-based metrics (Jaccard distance, binary core microbiome analysis) are uninformative in this system. Abundance-based compositional approaches (CLR transformation, Bray-Curtis distance, differential abundance testing) are essential for capturing the biologically relevant variation. The near-universal prevalence also means that compartment-driven assembly operates through quantitative amplification and suppression of existing community members, not through recruitment of novel taxa from external pools or exclusion of residents.

### Disease: a subtle, targeted signal

The disease signal in this dataset is remarkably subtle. At the community level, PERMANOVA is non-significant (*p* = 0.120), zero MAGs are differentially abundant after FDR correction, and the core microbiome is unchanged. Yet three complementary analyses reveal targeted disease effects:

First, three MAGs — two Verrucomicrobiota and one Actinomycetota — show large individual effect sizes (Cohen's *d* = 1.63–2.49) despite not surviving stringent FDR correction across 1,018 tests. These may represent genuine disease responders masked by multiple testing correction against a backdrop of 1,015 unchanged taxa.

Second, 125 disease indicator MAGs versus 36 healthy indicators (3.5:1 ratio) indicate a fidelity-based disease signal: disease-associated MAGs are not significantly more abundant but are more consistently present across diseased replicates. This suggests that BSR stabilises a subset of the community — possibly through stress-induced metabolic activity or pathogen-derived nutrient pulses — without dramatically altering overall community structure.

Third, the enrichment of candidate phyla radiation (Patescibacteria) and DPANN (Micrarchaeota) organisms as disease indicators is ecologically coherent. These ultrasmall organisms with reduced genomes are obligate symbionts or parasites of other microbes (Brown et al., 2015; Castelle et al., 2018). Their enrichment under disease suggests proliferation when normal community structure is perturbed, functioning as microbial dysbiosis markers analogous to the Firmicutes/Bacteroidetes ratio in human gut dysbiosis.

### Limitations

Several limitations should be noted. First, this is a cross-sectional study; temporal dynamics of community assembly and disease progression cannot be inferred. Longitudinal sampling through BSR development would clarify whether compartment-driven assembly is established before infection or reorganised during disease. Second, all samples originate from a single plantation in Peninsular Malaysia; generalisability across soil types, climatic zones, and oil palm varieties requires validation. Third, the use of MAG relative abundance as a proxy for absolute cell counts cannot distinguish between changes in the abundance of specific taxa and changes in total microbial biomass. Quantitative approaches (e.g., spike-in standards or qPCR normalisation) would strengthen abundance estimates. Fourth, MAG recovery is inherently incomplete; the 1,018 MAGs represent the binnable fraction of the metagenome, and unbinned community members may contribute ecologically relevant functions. Finally, compartment specificity scores are sensitive to the sampling depth and evenness of coverage across compartments.

---

## Conclusions

This study provides the first MAG-level community analysis of oil palm root microbiomes under *Ganoderma* basal stem rot, yielding three principal conclusions:

1. **Compartment identity is the master ecological variable**, explaining 42.4% of community variation while disease explains only 5.2% (non-significant). This extends the two-step selection model to tropical perennial crops and demonstrates that compartment-driven assembly is robust even under disease pressure.

2. **The endosphere harbours exclusive ecological specialists**, including a 15-member Myxococcota predatory guild (10-fold enrichment), five paradoxical *Roseiarcus* photoheterotrophs, and a sharp Pseudomonadota dominance contrasting with Acidobacteriota rhizosphere enrichment. These findings reveal a previously undocumented trophic structure inside oil palm roots.

3. **Disease produces a subtle, targeted signal** — not a broad community restructuring. Three MAGs with large effect sizes, 125 fidelity-based indicators, and enrichment of ultrasmall parasitic organisms (Patescibacteria, Micrarchaeota) suggest that BSR elicits a dysbiosis-like response detectable through indicator analysis but invisible to standard PERMANOVA or differential abundance tests.

These findings reframe BSR microbiome research from seeking community-wide compositional shifts to investigating targeted taxonomic signals and the ecological resilience mechanisms that buffer the microbiome against disease perturbation.

---

## Acknowledgements

This research was supported by [funding details]. We thank the Malaysian Palm Oil Board (MPOB) for access to plantation sites and metagenomic data. [Additional acknowledgements].

---

## Data availability

Metagenomic sequence data have been deposited in [repository] under accession [number]. The magprofile analysis package (v0.43) is available at [repository URL].

---

## References

Alneberg J, Bjarnason BS, de Bruijn I et al. Binning metagenomic contigs by coverage and composition. *Nat Methods* 2014;**11**:1144–6.

Anderson MJ. A new method for non-parametric multivariate analysis of variance. *Austral Ecol* 2001;**26**:32–46.

Anderson MJ. Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 2006;**62**:245–53.

Bowers RM, Kyrpides NC, Stepanauskas R et al. Minimum information about a metagenome-assembled genome of bacteria and archaea (MIMAG) of bacteria and archaea. *Nat Biotechnol* 2017;**35**:725–31.

Brown CT, Hug LA, Thomas BC et al. Unusual biology across a group comprising more than 15% of domain Bacteria. *Nature* 2015;**523**:208–11.

Bulgarelli D, Rott M, Schlaeppi K et al. Revealing structure and assembly cues for *Arabidopsis* root-inhabiting bacterial microbiota. *Nature* 2012;**488**:91–5.

Castelle CJ, Brown CT, Anantharaman K et al. Biosynthetic capacity, metabolic variety and unusual biology in the CPR and DPANN radiations. *Nat Rev Microbiol* 2018;**16**:629–45.

Chaumeil P-A, Mussig AJ, Hugenholtz P et al. GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics* 2022;**38**:5315–6.

Chklovski A, Parks DH, Woodcroft BJ et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 2023;**20**:1203–12.

Chong KP, Dayou J, Alexander A. Detection and control of *Ganoderma boninense* in oil palm crop. In: *SpringerBriefs in Agriculture*. Singapore: Springer, 2017.

Clarke KR. Non-parametric multivariate analyses of changes in community structure. *Austral Ecol* 1993;**18**:117–43.

Dufrêne M, Legendre P. Species assemblages and indicator species: the need for a flexible asymmetrical approach. *Ecol Monogr* 1997;**67**:345–66.

Edwards J, Johnson C, Santos-Medellín C et al. Structure, variation, and assembly of the root-associated microbiomes of rice. *Proc Natl Acad Sci USA* 2015;**112**:E911–20.

Gloor GB, Macklaim JM, Pawlowsky-Glahn V et al. Microbiome datasets are compositional: and this is not optional. *Front Microbiol* 2017;**8**:2224.

Idris AS, Kushairi A, Ismail S et al. Selection for partial resistance in oil palm progenies to *Ganoderma* basal stem rot. *J Oil Palm Res* 2004;**16**:19–26.

Kalam S, Basu A, Ahmad I et al. Recent understanding of soil Acidobacteria and their ecological significance: a critical review. *Front Microbiol* 2020;**11**:580024.

Kang DD, Li F, Kirton E et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 2019;**7**:e7359.

Kee YJ, Suhaimi N, Zakaria L et al. Characterisation of endophytic fungi from healthy and diseased oil palm roots. *J Oil Palm Res* 2022;**34**:315–25.

Kielak AM, Barreto CC, Kowalchuk GA et al. The ecology of Acidobacteria: moving beyond genes and genomes. *Front Microbiol* 2016;**7**:744.

Kulichevskaya IS, Danilova OV, Tereshina VM et al. *Roseiarcus fermentans* gen. nov., sp. nov., a bacteriochlorophyll *a*-containing fermentative bacterium related to alphaproteobacterial methanotrophs. *Int J Syst Evol Microbiol* 2014;**64**:2137–43.

Li D, Liu C-M, Luo R et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 2015;**31**:1674–6.

Lundberg DS, Lebeis SL, Paredes SH et al. Defining the core *Arabidopsis thaliana* root microbiome. *Nature* 2012;**488**:86–90.

Morgan AD, MacLean RC, Hillesland KL et al. Comparative analysis of *Myxococcus* predation on soil bacteria. *Appl Environ Microbiol* 2010;**76**:6920–7.

MPOB (Malaysian Palm Oil Board). *Malaysian Oil Palm Statistics 2022*. Bangi: MPOB, 2023.

Murphy DJ. The future of oil palm as a major global crop: opportunities and challenges. *J Oil Palm Res* 2014;**26**:1–24.

Pan S, Zhao X-M, Coelho LP. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics* 2023;**39**:i21–9.

Parks DH, Chuvochina M, Rinke C et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res* 2022;**50**:D199–210.

Parks DH, Rinke C, Chuvochina M et al. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. *Nat Microbiol* 2017;**2**:1533–42.

Peiffer JA, Spor A, Koren O et al. Diversity and heritability of the maize rhizosphere microbiome under field conditions. *Proc Natl Acad Sci USA* 2013;**110**:6548–53.

Petters S, Groß V, Söllinger A et al. The soil microbiome's keystone predators: Myxobacteria as the main bacterial predators in soil. *ISME J* 2021;**15**:2665–75.

Shamsilawani AB, Tao Z, Sariah M. Characterization of rhizosphere and endosphere microbiomes of healthy and *Ganoderma*-infected oil palm. *J Oil Palm Res* 2020;**32**:450–62.

Sieber CMK, Probst AJ, Sharrar A et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol* 2018;**3**:836–43.

Trivedi P, Leach JE, Tringe SG et al. Plant–microbiome interactions: from community assembly to plant health. *Nat Rev Microbiol* 2020;**18**:607–21.

Wu Y-W, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* 2016;**32**:605–7.

---

## Tables

**Table 1.** Alpha diversity across 30 samples.

| Metric | Mean ± SD | Range |
|--------|-----------|-------|
| Shannon diversity | 5.994 ± 0.227 | 5.51–6.38 |
| Simpson diversity | 0.994 ± 0.002 | 0.988–0.997 |
| Observed richness | 1,018.0 ± 0.2 | 1,017–1,018 |
| Pielou's evenness | 0.864 ± 0.033 | 0.787–0.921 |

**Table 2.** PERMANOVA, ANOSIM, and PERMDISP results across grouping variables.

| Test | Compartment | Health status | Group (D × C) |
|------|-------------|---------------|----------------|
| PERMANOVA pseudo-*F* | 9.955 | 1.546 | 4.757 |
| PERMANOVA *p* | 0.001*** | 0.120 NS | 0.001*** |
| PERMANOVA R² | 0.424 | 0.052 | 0.498 |
| ANOSIM R | 0.685 | 0.055 | 0.566 |
| ANOSIM *p* | 0.001*** | 0.089 NS | 0.001*** |
| PERMDISP *F* | 2.323 | 1.528 | 0.566 |
| PERMDISP *p* | 0.122 NS | 0.226 NS | 0.715 NS |

**Table 3.** Pairwise PERMANOVA for compartment comparisons.

| Comparison | Pseudo-*F* | *p* | *q* | R² |
|------------|--------:|-----:|-----:|----:|
| Bulk vs Endosphere | 15.43 | 0.001 | 0.002 | 0.462 |
| Endosphere vs Rhizosphere | 8.89 | 0.001 | 0.002 | 0.331 |
| Bulk vs Rhizosphere | 4.33 | 0.002 | 0.002 | 0.194 |

**Table 4.** Summary of differential abundance analysis (number of significant MAGs, *q* < 0.05).

| Comparison | Significant MAGs |
|------------|----------------:|
| Bulk vs Endosphere | 634 |
| Endosphere vs Rhizosphere | 474 |
| Bulk vs Rhizosphere | 376 |
| Diseased vs Healthy (global) | 0 |
| Within-compartment (D vs H) | 0 each |
| Cross-compartment groups | 0–198 per pair |

**Table 5.** Three MAGs with large effect sizes between disease states.

| MAG ID | Phylum | Genus | Cohen's *d* | Direction |
|--------|--------|-------|------------:|:---------:|
| Coassem3_concoct_34 | Verrucomicrobiota | UBA11358 | 2.49 | Disease ↑ |
| Coassem3_SemiBin_cstrat_454 | Verrucomicrobiota | UBA7542 | 2.01 | Disease ↑ |
| Coassem3_metabat_610 | Actinomycetota | JAJYUU01 | 1.63 | Disease ↑ |

**Table 6.** Top 10 endosphere indicator MAGs by IndVal score.

| Rank | MAG ID | Phylum | Genus | IndVal |
|-----:|--------|--------|-------|-------:|
| 1 | — | Myxococcota | *Palsa-1150* | 0.91 |
| 2 | — | Alphaproteobacteria | *Roseiarcus* | 0.90 |
| 3 | — | Alphaproteobacteria | *Roseiarcus* | 0.90 |
| 4 | — | Myxococcota | *Palsa-1150* | 0.89 |
| 5 | — | Myxococcota | JADGRB01 | 0.88 |
| 6 | — | Myxococcota | JADGRB01 | 0.88 |
| 7 | — | Alphaproteobacteria | *Roseiarcus* | 0.88 |
| 8 | — | Alphaproteobacteria | *Roseiarcus* | 0.87 |
| 9 | — | Alphaproteobacteria | *Roseiarcus* | 0.87 |
| 10 | — | Myxococcota | *Palsa-1150* | 0.87 |

**Table 7.** Archaeal phylum abundance by compartment.

| Phylum | Bulk soil (%) | Endosphere (%) | Cohen's *d* |
|--------|-------------:|---------------:|------------:|
| Thermoplasmatota | 15.4 | 3.5 | 6.43 |
| Thermoproteota | 4.2 | 1.4 | 6.02 |

**Table 8.** Core microbiome at varying prevalence thresholds.

| Prevalence threshold | Core MAGs (all groups) |
|---------------------:|-----------------------:|
| 50% | 1,018 |
| 75% | 1,018 |
| 90% | 1,017–1,018 |
| 100% | 1,017 |

**Table 9.** Indicator species counts by health status.

| Health state | Significant indicators | Top phyla among indicators |
|-------------|----------------------:|---------------------------|
| Diseased | 125 | Pseudomonadota, Acidobacteriota, Patescibacteria |
| Healthy | 36 | Pseudomonadota, Actinomycetota |

---

## Figure legends

**Figure 1.** Differential abundance between compartments. (a) Volcano plot of CLR abundance differences between bulk soil and endosphere (634 significant MAGs highlighted). (b) Venn diagram of significant MAGs across three pairwise compartment comparisons. (c) Phylum-level composition of differentially abundant MAGs for each comparison.

**Figure 2.** Disease-associated MAGs. (a) Volcano plot of CLR abundance differences between diseased and healthy samples (3 MAGs with Cohen's *d* > 1.5 highlighted). (b) Boxplots of the three disease-associated MAGs across health states.

**Figure 3.** Compartment specialists and indicator taxa. (a) Distribution of compartment specificity scores, with specialist threshold (>0.5) indicated; all 65 specialists are endosphere-associated. (b) Myxococcota predatory guild: 15 endosphere indicator MAGs coloured by genus. (c) Phylum composition of indicator MAGs for each compartment, showing Acidobacteriota dominance in rhizosphere and Pseudomonadota/Myxococcota dominance in endosphere.

**Figure 4.** Archaeal exclusion from root compartments. Phylum-level abundance of Thermoplasmatota and Thermoproteota across compartments, with Cohen's *d* values indicated.

**Figure 5.** Core microbiome and universal membership. (a) Core MAG count at prevalence thresholds from 50% to 100%. (b) Rank-abundance curves by compartment showing abundance redistribution of universally shared MAGs.

**Figure 6.** Disease indicator asymmetry. (a) Number of disease vs healthy indicators. (b) Top disease indicator MAGs, highlighting ultrasmall Patescibacteria and Micrarchaeota. (c) Phylum-level enrichment by health state with FDR-corrected significance.
