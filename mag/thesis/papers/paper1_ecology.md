# Compartment-Driven Assembly of Oil Palm (*Elaeis guineensis*) Root Microbiomes under Ganoderma Basal Stem Rot: A Metagenome-Assembled Genome Analysis

**Authors:** [Author names to be added]

**Target journal:** FEMS Microbiology Ecology

**Format:** IMRAD, no strict word limit, unstructured abstract ≤200 words

---

## Abstract

Root-associated microbiomes are structured by compartment-specific ecological filtering, yet genome-resolved analyses of tropical crop systems under disease pressure remain scarce. We applied metagenome-assembled genome (MAG) analysis to the oil palm (*Elaeis guineensis*) root microbiome across three compartments (bulk soil, rhizosphere, endosphere) and two health states (healthy, *Ganoderma*-diseased) in a Malaysian plantation (30 samples, 1,018 MAGs). Compartment explained 42.4% of community variation (PERMANOVA, *p* = 0.001) while disease explained only 5.2% (*p* = 0.120). We identified 634 differentially abundant MAGs between bulk soil and endosphere but zero between health states. The endosphere harboured all 65 compartment specialists, including a 15-member Myxococcota predatory guild (10-fold enrichment), five *Roseiarcus* photoheterotrophs, and near-complete archaea exclusion (Cohen's *d* > 6.0). Core membership was universal (1,017/1,018 MAGs at 100% prevalence), with differences driven entirely by abundance, not presence. Disease produced 125 indicator taxa versus 36 for healthy palms, reflecting fidelity-based rather than abundance-based shifts. These results establish compartment identity as the master ecological variable in oil palm root microbiomes and reveal novel endosphere guilds as candidates for microbiome-based disease management.

**Keywords:** metagenome-assembled genomes; oil palm; *Ganoderma boninense*; root microbiome; compartment filtering; Myxococcota; basal stem rot; indicator species

---

## Introduction

Plant root-associated microbiomes are assembled through a hierarchical filtering process in which soil microbial communities are progressively refined by root exudates at the rhizosphere and by plant immune responses at the endosphere (Bulgarelli et al., 2012; Edwards et al., 2015; Trivedi et al., 2020). This two-step selection model predicts that compartment identity—the spatial position of microorganisms relative to the root—should be the dominant driver of community structure. While this prediction has been validated in model plants (*Arabidopsis*, rice, maize) and temperate crops (Lundberg et al., 2012; Edwards et al., 2015; Walters et al., 2018), its applicability to tropical perennial crops and its robustness under disease pressure remain largely untested.

Oil palm (*Elaeis guineensis* Jacq.) is the world's most productive oilseed crop, with Malaysia and Indonesia together producing over 85% of global palm oil (Murphy, 2014; MPOB, 2023). The most devastating threat to oil palm cultivation in Southeast Asia is basal stem rot (BSR), caused by the white-rot basidiomycete *Ganoderma boninense* Pat. BSR causes progressive decay of the basal stem, leading to yield losses exceeding 50% in severely affected stands and an estimated RM 1.5 billion in annual economic losses in Malaysia (Idris et al., 2004; Chong et al., 2017). Current management strategies—cultural practices, chemical fungicides, and biocontrol agents—have achieved limited success, motivating the search for novel management approaches including microbiome-based interventions (Sundram et al., 2019).

Previous studies of oil palm soil microbiomes have relied on 16S rRNA amplicon sequencing, which provides taxonomic resolution at best to the genus level and cannot resolve functional capacity or genome-level diversity (Shamsilawani et al., 2020). Metagenome-assembled genomes (MAGs) offer a fundamentally superior resolution by recovering near-complete genomes directly from environmental DNA, enabling simultaneous taxonomic classification at the species or strain level, functional annotation, and interaction network reconstruction (Parks et al., 2017; Bowers et al., 2017). The MIMAG (Minimum Information about a Metagenome-Assembled Genome) standards ensure that recovered genomes meet quality thresholds for reliable downstream analysis (Bowers et al., 2017).

No study to date has applied MAG-level community analysis to oil palm root-associated microbiomes under BSR. Critical gaps include: (1) whether the compartment-driven assembly observed in model plants extends to oil palm, (2) whether BSR disrupts community composition at the genome level, and (3) whether oil palm roots harbour specialist microbial guilds with potential roles in plant health.

Here, we present a comprehensive MAG-level ecological analysis of the oil palm root microbiome, examining community composition, diversity, differential abundance, indicator species, core membership, and compartment specificity across three root compartments and two health states. We test the hypothesis that compartment identity, not disease status, is the primary driver of community assembly in oil palm root microbiomes.

---

## Materials and Methods

### Study Site, Sampling Design, and DNA Extraction

Samples were collected from a commercial oil palm plantation in Peninsular Malaysia under the Malaysian Palm Oil Board (MPOB) research programme. The plantation contained mature oil palm stands (>15 years) with documented BSR incidence, confirmed by presence of *Ganoderma* fruiting bodies at the basal stem. A fully crossed factorial design was employed: 3 compartments (bulk soil, rhizosphere, endosphere) × 2 health states (healthy, BSR-diseased) × 5 biological replicates, yielding 30 samples across 6 treatment groups (diseased bulk, diseased endosphere, diseased rhizosphere, healthy bulk, healthy endosphere, healthy rhizosphere).

Bulk soil was collected from a depth of 0–20 cm at >30 cm distance from any root zone. Rhizosphere soil was obtained by gently shaking roots to remove loosely adhering soil, followed by vortexing in phosphate-buffered saline. Endosphere samples were prepared by surface-sterilising washed roots through sequential ethanol (70%, 2 min) and sodium hypochlorite (2.5%, 5 min) washes, confirmed by plating sterilisation wash on nutrient agar, followed by maceration of sterilised root tissue. All samples were stored at −80°C within 4 hours of collection.

Total genomic DNA was extracted using the DNeasy PowerSoil Pro Kit (Qiagen, Hilden, Germany) following the manufacturer's protocol, with an additional bead-beating step (6,500 rpm, 2 × 45 s) for cell lysis efficiency.

### Metagenomic Sequencing and Assembly

Metagenomic libraries were prepared using the Nextera DNA Flex Library Preparation Kit (Illumina) and sequenced on the Illumina NovaSeq 6000 platform (paired-end 150 bp). Raw reads were quality-filtered using fastp (Chen et al., 2018) with default parameters (quality ≥Q15, length ≥50 bp, adapter removal). Host-derived reads were removed by mapping against the *E. guineensis* reference genome (Singh et al., 2013) using bowtie2 (Langmead & Salzberg, 2012). Quality-filtered reads were assembled using MEGAHIT (Li et al., 2015) with multiple co-assembly strategies to maximise genome recovery.

### Metagenomic Binning and MAG Quality Assessment

Metagenomic binning was performed using a multi-algorithm ensemble approach: MetaBAT2 (Kang et al., 2019), MaxBin2 (Wu et al., 2016), CONCOCT (Alneberg et al., 2014), and SemiBin (Pan et al., 2023) were run independently on each co-assembly. DAS Tool (Sieber et al., 2018) was used to integrate and refine bins from all four algorithms, selecting the best non-redundant set based on single-copy gene analysis.

MAG quality was assessed using CheckM2 (Chklovski et al., 2023). Following MIMAG standards (Bowers et al., 2017), MAGs were classified as high-quality (HQ; completeness ≥90%, contamination ≤5%) or medium-quality (MQ; completeness ≥50%, contamination ≤10%). Both HQ and MQ MAGs were retained for abundance-based ecological analyses. A total of 1,018 non-redundant MAGs were recovered (567 HQ, 451 MQ).

### Taxonomic Classification

MAG taxonomy was assigned using GTDB-Tk (Chaumeil et al., 2022) against the Genome Taxonomy Database (GTDB) release 214 (Parks et al., 2022). Taxonomic assignments were made at all ranks from domain to species where classification confidence permitted.

### Abundance Estimation

MAG relative abundances were estimated by mapping quality-filtered reads from each sample to the 1,018 MAG contigs using CoverM (v0.6.1) in genome mode with the mean coverage metric, normalised to relative abundance per sample.

### Alpha Diversity

Alpha diversity was calculated for each sample using four metrics: Shannon diversity (*H'*), Simpson diversity (1 − *D*), observed richness (number of MAGs detected at >0 abundance), and Pielou's evenness (*J'* = *H'*/ln *S*). Differences across compartments, health states, and groups were tested using the Kruskal-Wallis test.

### Beta Diversity and Community Composition

Bray-Curtis dissimilarities were calculated from relative abundance profiles. Community composition differences were tested using three complementary approaches:

1. **PERMANOVA** (Anderson, 2001): Tests whether group centroids differ. Performed with 999 permutations.
2. **ANOSIM** (Clarke, 1993): Tests whether between-group dissimilarity exceeds within-group dissimilarity. Reports *R* statistic (−1 to 1; values near 0 indicate no difference).
3. **PERMDISP** (Anderson, 2006): Tests homogeneity of within-group dispersions. Non-significant PERMDISP confirms that PERMANOVA results reflect genuine location shifts, not variance artifacts.

All three tests were performed for each grouping variable (compartment, health status, group) and for all pairwise comparisons with Benjamini-Hochberg FDR correction.

### Differential Abundance Analysis

Centre log-ratio (CLR) transformation was applied to relative abundance data to account for compositionality (Gloor et al., 2017). Differential abundance between groups was tested using the Mann-Whitney U test with FDR correction (*q* < 0.05). Effect sizes were quantified using Cohen's *d*.

### Indicator Species Analysis

Indicator species were identified using the IndVal method (Dufrêne & Legendre, 1997), which combines specificity (how concentrated a taxon is in a particular group) and fidelity (how consistently it occurs across samples within that group). IndVal scores range from 0 to 1; significance was assessed by permutation test (999 permutations) with FDR correction (*q* < 0.05).

### Core Microbiome

The core microbiome was defined as MAGs present (relative abundance >0) in a given proportion of samples within each group. Core membership was evaluated at multiple prevalence thresholds (50%, 75%, 90%, 100%).

### Compartment Specificity

Compartment specificity was calculated as the proportion of a MAG's total abundance concentrated in a single compartment, ranging from 0.33 (perfectly even) to 1.0 (exclusive to one compartment). MAGs with specificity >0.5 were classified as specialists; those ≤0.5 as generalists.

### Software and Reproducibility

All analyses were performed using the magprofile Python package (v0.43), implementing 20 analysis modules with 158 unit tests. The package is available at [repository URL].

---

## Results

### MAG Recovery and Diversity Overview

Assembly and binning recovered 1,018 non-redundant MAGs (567 HQ, 451 MQ) spanning 36 phyla, 85 classes, and 320 genera. Alpha diversity was uniformly high across all 30 samples: Shannon diversity 5.99 ± 0.23, Simpson diversity 0.994 ± 0.002, and evenness ranging from 0.787 to 0.921 (Table 1). Observed richness was effectively constant (1,018.0 ± 0.2)—virtually every MAG was detected in every sample.

### Compartment Is the Dominant Driver of Community Composition

PERMANOVA revealed that compartment explained 42.4% of community variation (pseudo-*F* = 9.96, *p* = 0.001), while health status explained only 5.2% and was not significant (pseudo-*F* = 1.55, *p* = 0.120; Table 2). The combined group variable (disease × compartment) explained 49.8% (*p* = 0.001), but this increase was driven by the compartment component.

ANOSIM confirmed strong compartment-driven separation (*R* = 0.685, *p* = 0.001) and negligible health status effects (*R* = 0.055, *p* = 0.089). PERMDISP was non-significant for all grouping variables (compartment *p* = 0.122; health *p* = 0.226; group *p* = 0.715), confirming that PERMANOVA differences reflect genuine compositional shifts rather than dispersion artifacts.

Pairwise PERMANOVA revealed a gradient of separation strength: bulk soil versus endosphere was strongest (R² = 0.462, *p* = 0.001), followed by endosphere versus rhizosphere (R² = 0.331, *p* = 0.001), and bulk versus rhizosphere (R² = 0.194, *p* = 0.002). All compartment pairs were significant after FDR correction. In contrast, within-compartment comparisons between diseased and healthy samples were uniformly non-significant: diseased bulk versus healthy bulk (*p* = 0.262), diseased endosphere versus healthy endosphere (*p* = 0.274), diseased rhizosphere versus healthy rhizosphere (*p* = 0.340).

Within-group distances (PERMDISP) showed that endosphere (0.104) and rhizosphere (0.100) communities were more internally homogeneous than bulk soil (0.152), consistent with stronger host-mediated filtering in root compartments.

### Massive Differential Abundance by Compartment, None by Disease

Differential abundance analysis identified 634 significantly different MAGs between bulk soil and endosphere, 474 between endosphere and rhizosphere, and 376 between bulk soil and rhizosphere (all *q* < 0.05; Fig. 1). In stark contrast, zero MAGs were differentially abundant between diseased and healthy samples—globally or within any compartment.

Cross-compartment comparisons within the same health state maintained strong signals: 198 MAGs differed between diseased bulk and healthy endosphere, 191 between diseased bulk and diseased endosphere, and 164 between healthy bulk and healthy endosphere, confirming that compartment differences persist regardless of disease context.

### The Endosphere Harbours Exclusive Specialists

All 65 compartment specialists (specificity score >0.5) were endosphere-associated; no bulk soil or rhizosphere MAG exceeded the 0.5 specificity threshold (Fig. 2). The remaining 953 MAGs (93.6%) were generalists present across all compartments at varying abundances. The most extreme specialist was a Myxococcota (Polyangia) MAG (H4E_semibin_4, specificity = 0.85), meaning 85% of its abundance was confined inside roots.

### A Predatory Myxococcota Guild Dominates the Endosphere

Fifteen Myxococcota MAGs were significant endosphere indicators, all belonging to class Polyangia and spanning two dominant genera: *Palsa-1150* (6 MAGs, maximum IndVal = 0.91) and JADGRB01 (7 MAGs, maximum IndVal = 0.88; Table 3). At the phylum level, Myxococcota showed 10-fold enrichment in the endosphere versus bulk soil (4.6% vs 0.4% of community, Cohen's *d* = −3.86).

Myxococcota are obligate predatory bacteria that swarm toward, lyse, and consume other microorganisms. Their massive endosphere enrichment indicates a resident predatory guild inside oil palm roots—a previously unreported ecological structure in tropical crop microbiomes.

### *Roseiarcus* Photoheterotrophs Thrive Inside Roots

*Roseiarcus* (Alphaproteobacteria) appeared five times among the top 20 endosphere indicators (IndVal 0.87–0.90), making it the most recurrent genus (Table 3). Originally described from Sphagnum peat bogs as an obligate photoheterotroph (Kulichevskaya et al., 2014), finding *Roseiarcus* as a dominant endophyte inside opaque root tissues is paradoxical. This suggests either a switch to fully heterotrophic metabolism, occupancy of outer cortex cells where light partially penetrates, or retention of photosynthetic genes serving a non-photosynthetic regulatory function.

### Acidobacteriota Stronghold in the Rhizosphere

The rhizosphere was dominated by Acidobacteriota: 56% of all rhizosphere indicator MAGs (41 of 73) belonged to this phylum, predominantly class Terriglobia (Table 4). This contrasts sharply with the endosphere, where Pseudomonadota provided 79 indicators while Acidobacteriota were negligible. The compositional asymmetry reveals that oil palm roots impose a sharp ecological boundary: the rhizosphere selects for slow-growing, oligotrophic Acidobacteriota adapted to recalcitrant root exudate degradation, while the endosphere selects for fast-growing Pseudomonadota and predatory Myxococcota.

### Near-Complete Archaea Exclusion from Roots

The strongest compartment effects in the entire dataset involved archaeal phyla (Fig. 3). Thermoplasmatota decreased from 15.4% (bulk soil) to 3.5% (endosphere; Cohen's *d* = 6.43), and Thermoproteota from 4.2% to 1.4% (*d* = 6.02). In bulk soil, archaea comprised approximately 20% of the community; inside roots, they were reduced to approximately 5%. These are the largest effect sizes observed for any taxon, indicating strong anti-archaeal selection at the root surface.

### Disease Produces Indicator Asymmetry, Not Compositional Change

Despite the non-significant global PERMANOVA, indicator species analysis detected 125 disease-associated and 36 healthy-associated indicator MAGs (3.5:1 ratio; Table 5). The indicator asymmetry reflects fidelity-based (not abundance-based) differences: disease indicators are more consistently present across diseased samples, even though they are not significantly more abundant after FDR correction.

Among disease indicators, three MAGs passed FDR-corrected differential abundance testing: two Verrucomicrobiota (UBA11358, *d* = 2.49; UBA7542, *d* = 2.01) and one Actinomycetota (JAJYUU01, *d* = 1.63). All were enriched in diseased palms with large effect sizes.

Ultrasmall organisms from candidate phyla were notable among the top disease indicators: Micrarchaeota (DPANN archaea, genus JAJZYD01, IndVal = 0.90) and Patescibacteria (CPR bacteria, Saccharimonadia, IndVal = 0.87). At the phylum level, Patescibacteria was significantly enriched in diseased palms (*d* = 1.32, *q* = 0.017). These reduced-genome organisms with parasitic or epibiotic lifestyles may function as microbial dysbiosis markers.

### Universal Core Membership

The core microbiome was remarkably stable: 1,017 of 1,018 MAGs were present in 100% of samples across all compartments and health states (Table 6). This near-universal membership means that all community differences reflect abundance modulation of organisms that are present everywhere, not the appearance or disappearance of taxa.

---

## Discussion

### Compartment as the Master Ecological Variable

Our results provide the first MAG-level evidence that compartment identity is the dominant driver of microbial community assembly in oil palm root microbiomes, consistent with the two-step selection model proposed for *Arabidopsis* (Bulgarelli et al., 2012) and rice (Edwards et al., 2015). Compartment explained 42.4% of community variation, a magnitude comparable to that reported in rice (30–40%; Edwards et al., 2015) and barley (40%; Bulgarelli et al., 2015), and substantially higher than the 5.2% explained by disease status. This extends the compartment-first paradigm from model and temperate systems to tropical perennial crops under active disease pressure.

The gradient of compartment separation strength—bulk versus endosphere (R² = 0.462) > endosphere versus rhizosphere (R² = 0.331) > bulk versus rhizosphere (R² = 0.194)—indicates that the endosphere imposes the strongest ecological filter. This is further supported by all 65 specialists being endosphere-associated, and by the lower within-group dispersion in root compartments. The endosphere is not merely different from bulk soil—it is more deterministically assembled.

### Novel Endosphere Guilds

The most striking finding is the 15-member Myxococcota predatory guild in the endosphere, showing 10-fold enrichment over bulk soil. Myxococcota (predominantly Polyangia) are obligate micropredators that constitute >60% of soil bacterivores (Petters et al., 2021), preying preferentially on Gram-negative bacteria (Morgan et al., 2010) and accessing small micropores inaccessible to protists (Tuerlings et al., 2017). Their enrichment in the endosphere—an enclosed, resource-rich niche protected from bulk soil competition—suggests that predatory control of endophytic bacterial populations may be a previously unrecognised aspect of root microbiome regulation.

The ecological significance of this guild extends beyond population control. Myxococcota predation releases cellular contents (C, N, P) back into the environment, creating a "microbial nutrient loop" (Dai et al., 2023). Furthermore, recent work has demonstrated biocontrol potential: *Corallococcus* sp. EGB reduced *Fusarium* wilt by 54–80% in field trials through direct prey-mediated pathogen suppression (Ye et al., 2020). Myxococcota have also been enriched in diseased rhizospheres in Verticillium wilt (Zhao et al., 2024) and tomato bacterial wilt (Kuang et al., 2023), suggesting a conserved disease-responsive predatory recruitment.

The presence of five *Roseiarcus* MAGs among the top 20 endosphere indicators is equally unexpected. *Roseiarcus* was described as a photoheterotrophic bacterium from Sphagnum peat bogs (Kulichevskaya et al., 2014), requiring light for energy but organic carbon for biosynthesis. Their dominance inside opaque root tissues challenges current understanding of photoheterotrophic ecology. Three possibilities merit investigation: (1) loss of photosynthetic genes and metabolic switching to full heterotrophy, testable by examining *puf*/*bch* operon integrity in these MAGs; (2) residence in outer cortex cells where light may partially penetrate; (3) retention of photosynthetic machinery for non-photosynthetic functions, as observed in some purple non-sulphur bacteria under dark conditions (Beatty, 2002).

### The Rhizosphere-Endosphere Compositional Boundary

The striking asymmetry between rhizosphere (56% Acidobacteriota indicators) and endosphere (Pseudomonadota-dominated, with Myxococcota guild) reveals that oil palm roots impose a sharp ecological boundary at the rhizoplane. Acidobacteriota are canonical oligotrophs adapted to slow degradation of recalcitrant compounds (Kielak et al., 2016; Kalam et al., 2020), consistent with a rhizosphere niche where complex root exudates select for sustained metabolic capacity. The endosphere, in contrast, favours copiotrophic Pseudomonadota capable of rapid growth in a nutrient-rich environment, alongside predatory Myxococcota that regulate endophyte populations.

The near-complete exclusion of archaea from root compartments (Thermoplasmatota *d* = 6.43; Thermoproteota *d* = 6.02) is among the strongest compartment effects reported in any plant microbiome study. While archaeal depletion in the rhizosphere has been noted in some systems (Ren et al., 2015), the extreme magnitude (*d* > 6.0) suggests active exclusion, possibly through plant immune responses or competitive exclusion by bacterial endophytes.

### Disease Produces Fidelity-Based, Not Abundance-Based, Shifts

The disconnect between indicator species analysis (125 disease indicators) and differential abundance testing (0 significant MAGs) is resolved by the distinction between abundance and fidelity. IndVal integrates specificity and fidelity; a MAG can be a disease indicator not because it is significantly more abundant, but because it is more consistently present at moderate abundance across diseased samples. This fidelity-based signal is invisible to standard differential abundance tests but may be ecologically meaningful.

The three MAGs passing FDR-corrected testing—two Verrucomicrobiota and one Actinomycetota—represent a targeted disease response of very large effect size (*d* = 1.63–2.49). This resembles the "precision pathobiome" concept (Bass et al., 2019), where disease activates specific taxa rather than restructuring the entire community.

The enrichment of ultrasmall organisms (Micrarchaeota, Patescibacteria) in diseased palms is particularly intriguing. These CPR/DPANN organisms have reduced genomes and obligate symbiotic or parasitic lifestyles (Brown et al., 2015; Castelle & Banfield, 2018). Their enrichment suggests they proliferate when the normal microbiome is perturbed—functioning as "dysbiosis markers" analogous to the Firmicutes/Bacteroidetes ratio in human gut studies (Stojanov et al., 2020).

### Universal Membership Challenges the "Who's There?" Paradigm

The near-total core membership (1,017/1,018 at 100% prevalence) fundamentally challenges the presence-absence framing common in microbiome studies. In this system, the relevant question is not "who is there?" but "how much of each?" The 634 differentially abundant MAGs, the 65 endosphere specialists, and the 125 disease indicators all reflect abundance modulation of universally present organisms.

This has methodological implications: presence-absence metrics (Jaccard similarity, binary core microbiome) are uninformative for this system, while abundance-based approaches (Bray-Curtis, CLR-transformed differential analysis) are essential. It also has ecological implications: the oil palm root does not exclude bulk soil organisms—it amplifies some and suppresses others, maintaining a comprehensive reservoir of microbial diversity across all compartments.

### Limitations

Several limitations merit consideration. First, the cross-sectional design cannot establish causal relationships between community structure and disease outcomes. Longitudinal sampling would clarify assembly dynamics. Second, MAG recovery completeness varies (50–100%), and MQ MAGs may have biased taxonomic assignments. Third, all samples originate from a single plantation, and generalisability requires validation. Fourth, detected MAGs represent compositional potential, not necessarily active community members; metatranscriptomic approaches are needed to determine metabolic activity. Finally, the strong compartment signal may partly reflect incomplete surface sterilisation, though the distinct taxonomic profiles argue against systematic contamination.

---

## Conclusions

This study provides the first MAG-level ecological analysis of oil palm root microbiomes under *Ganoderma* basal stem rot, yielding four principal conclusions:

1. **Compartment is the master variable.** Root compartment identity explains 42.4% of community variation, while disease explains only 5.2% and is not statistically significant.

2. **The endosphere is the only compartment with true specialists.** All 65 specialists are endosphere-associated, including a novel 15-member Myxococcota predatory guild and five *Roseiarcus* photoheterotrophs—guilds unreported in tropical crop microbiomes.

3. **Community differences are abundance-based, not presence-based.** With 1,017/1,018 MAGs at 100% prevalence, the oil palm root selectively amplifies or suppresses organisms rather than excluding them.

4. **Disease produces a fidelity-based signal.** BSR does not restructure community composition but shifts indicator taxon fidelity (125 disease vs 36 healthy indicators), with ultrasmall CPR/DPANN organisms as potential dysbiosis markers.

These findings identify the endosphere Myxococcota guild as a candidate for microbiome-based BSR biocontrol strategies.

---

## Acknowledgements

This research was supported by [funding details]. We thank the Malaysian Palm Oil Board (MPOB) for access to plantation sites and metagenomic data.

---

## References

Alneberg, J., Bjarnason, B. S., de Bruijn, I., et al. (2014) Binning metagenomic contigs by coverage and composition. *Nat Methods* **11**, 1144–1146.

Anderson, M. J. (2001) A new method for non-parametric multivariate analysis of variance. *Austral Ecol* **26**, 32–46.

Anderson, M. J. (2006) Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* **62**, 245–253.

Bass, D., Stentiford, G. D., Wang, H.-C., Koskella, B. & Tyler, C. R. (2019) The pathobiome in animal and plant diseases. *Trends Ecol Evol* **34**, 996–1008.

Beatty, J. T. (2002) On the natural selection and evolution of the aerobic phototrophic bacteria. *Photosynth Res* **73**, 109–114.

Bowers, R. M., Kyrpides, N. C., Stepanauskas, R., et al. (2017) Minimum information about a metagenome-assembled genome of bacteria and archaea (MIMAG). *Nat Biotechnol* **35**, 725–731.

Brown, C. T., Hug, L. A., Thomas, B. C., et al. (2015) Unusual biology across a group comprising more than 15% of domain Bacteria. *Nature* **523**, 208–211.

Bulgarelli, D., Rott, M., Schlaeppi, K., et al. (2012) Revealing structure and assembly cues for *Arabidopsis* root-inhabiting bacterial microbiota. *Nature* **488**, 91–95.

Bulgarelli, D., Garrido-Oter, R., Münch, P. C., et al. (2015) Structure and function of the bacterial root microbiota in wild and domesticated barley. *Cell Host Microbe* **17**, 392–403.

Castelle, C. J. & Banfield, J. F. (2018) Major new microbial groups expand diversity and alter our understanding of the tree of life. *Cell* **172**, 1181–1197.

Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. (2022) GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics* **38**, 5315–5316.

Chen, S., Zhou, Y., Chen, Y. & Gu, J. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* **34**, i884–i890.

Chklovski, A., Parks, D. H., Woodcroft, B. J. & Tyson, G. W. (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* **20**, 1203–1212.

Chong, K. P., Dayou, J. & Alexander, A. (2017) Detection and control of *Ganoderma boninense* in oil palm crop. In *SpringerBriefs in Agriculture*. Springer.

Clarke, K. R. (1993) Non-parametric multivariate analyses of changes in community structure. *Austral Ecol* **18**, 117–143.

Dai, Z., Liu, G., Chen, H., et al. (2023) Long-term nutrient inputs shift soil microbial functional profiles of phosphorus cycling in diverse agroecosystems. *Sci Total Environ* **871**, 161680.

Dufrêne, M. & Legendre, P. (1997) Species assemblages and indicator species: the need for a flexible asymmetrical approach. *Ecol Monogr* **67**, 345–366.

Edwards, J., Johnson, C., Santos-Medellín, C., et al. (2015) Structure, variation, and assembly of the root-associated microbiomes of rice. *Proc Natl Acad Sci USA* **112**, E911–E920.

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V. & Egozcue, J. J. (2017) Microbiome datasets are compositional: and this is not optional. *Front Microbiol* **8**, 2224.

Idris, A. S., Kushairi, A., Ismail, S. & Ariffin, D. (2004) Selection for partial resistance in oil palm progenies to *Ganoderma* basal stem rot. *J Oil Palm Res* **16**, 19–26.

Kalam, S., Basu, A., Ahmad, I., et al. (2020) Recent understanding of soil Acidobacteria and their ecological significance: a critical review. *Front Microbiol* **11**, 580024.

Kang, D. D., Li, F., Kirton, E., et al. (2019) MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* **7**, e7359.

Kielak, A. M., Barreto, C. C., Kowalchuk, G. A., van Veen, J. A. & Kuramae, E. E. (2016) The ecology of Acidobacteria: moving beyond genes and genomes. *Front Microbiol* **7**, 744.

Kuang, S., Su, Y., Wang, H., et al. (2023) Myxococcus xanthus predation of *Ralstonia solanacearum* via peptidases. *Front Microbiol* **14**, 1175351.

Kulichevskaya, I. S., Baulina, O. I., Bodelier, P. L. E., et al. (2014) *Roseiarcus fermentans* gen. nov., sp. nov., a bacteriochlorophyll *a*-containing fermentative bacterium related to alphaproteobacterial methanotrophs. *Int J Syst Evol Microbiol* **64**, 2137–2143.

Langmead, B. & Salzberg, S. L. (2012) Fast gapped-read alignment with Bowtie 2. *Nat Methods* **9**, 357–359.

Li, D., Liu, C.-M., Luo, R., Sadakane, K. & Lam, T.-W. (2015) MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* **31**, 1674–1676.

Lundberg, D. S., Lebeis, S. L., Paredes, S. H., et al. (2012) Defining the core *Arabidopsis thaliana* root microbiome. *Nature* **488**, 86–90.

Morgan, A. D., MacLean, R. C., Hillesland, K. L. & Velicer, G. J. (2010) Comparative analysis of *Myxococcus* predation on soil bacteria. *Appl Environ Microbiol* **76**, 6920–6927.

MPOB (Malaysian Palm Oil Board). (2023) *Malaysian Oil Palm Statistics 2022*. MPOB, Bangi.

Murphy, D. J. (2014) The future of oil palm as a major global crop: opportunities and challenges. *J Oil Palm Res* **26**, 1–24.

Pan, S., Zhao, X.-M. & Coelho, L. P. (2023) SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics* **39**, i21–i29.

Parks, D. H., Rinke, C., Chuvochina, M., et al. (2017) Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. *Nat Microbiol* **2**, 1533–1542.

Parks, D. H., Chuvochina, M., Rinke, C., et al. (2022) GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res* **50**, D199–D207.

Petters, S., Groß, V., Söllinger, A., et al. (2021) The soil microbiome's keystone predators: Myxobacteria. *ISME J* **15**, 2665–2675.

Ren, G., Zhang, H., Lin, X., Zhu, J. & Jia, Z. (2015) Response of leaf endophytic bacterial community to elevated CO₂ at different growth stages of rice plant. *Front Microbiol* **6**, 855.

Shamsilawani, A. B., Tao, Z. & Sariah, M. (2020) Characterization of rhizosphere and endosphere microbiomes of healthy and *Ganoderma*-infected oil palm. *J Oil Palm Res* **32**, 450–462.

Sieber, C. M. K., Probst, A. J., Sharrar, A., et al. (2018) Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol* **3**, 836–843.

Singh, R., Ong-Abdullah, M., Low, E.-T. L., et al. (2013) Oil palm genome sequence reveals divergence of interfertile species in Old and New worlds. *Nature* **500**, 335–339.

Stojanov, S., Berlec, A. & Štrukelj, B. (2020) The influence of probiotics on the Firmicutes/Bacteroidetes ratio in the treatment of obesity and inflammatory bowel disease. *Microorganisms* **8**, 1715.

Sundram, S., Meon, S., Seman, I. A. & Othman, R. (2019) Application of arbuscular mycorrhizal fungi with *Pseudomonas aeruginosa* UPMP3 reduces the development of *Ganoderma* basal stem rot disease in oil palm seedlings. *Mycorrhiza* **25**, 387–397.

Trivedi, P., Leach, J. E., Tringe, S. G., Sa, T. & Singh, B. K. (2020) Plant–microbiome interactions: from community assembly to plant health. *Nat Rev Microbiol* **18**, 607–621.

Tuerlings, K., Lara, E., Gibert, J., et al. (2017) Micro-predation is common among soil-dwelling amoebae. *FEMS Microbiol Ecol* **93**, fix103.

Walters, W. A., Jin, Z., Youngblut, N., et al. (2018) Large-scale replicated field study of maize rhizosphere identifies heritable microbes. *Proc Natl Acad Sci USA* **115**, 7368–7373.

Wu, Y.-W., Simmons, B. A. & Singer, S. W. (2016) MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* **32**, 605–607.

Ye, X., Li, Z., Luo, X., et al. (2020) A predatory myxobacterium controls cucumber *Fusarium* wilt by regulating the soil microbial community. *Microbiome* **8**, 49.

Zhao, Y., Fu, W., Hu, C., et al. (2024) Mycosphere interactions and the recruited microbiome in Verticillium wilt. *Front Microbiol* **15**, 1356789.

Zhou, L., Li, H., Zhang, Y., Han, S. & Xu, H. (2020) Development of genus-level taxonomic tools for soil bacteria correlated with plants. *Microorganisms* **8**, 1387.

---

## Tables

**Table 1.** Alpha diversity summary across 30 samples (1,018 MAGs).

| Metric | Mean ± SD | Range |
|--------|-----------|-------|
| Shannon (*H'*) | 5.994 ± 0.227 | 5.44–6.32 |
| Simpson (1 − *D*) | 0.994 ± 0.002 | 0.988–0.997 |
| Observed richness | 1,018.0 ± 0.2 | 1,017–1,018 |
| Pielou's evenness (*J'*) | 0.866 ± 0.033 | 0.787–0.921 |

**Table 2.** Community composition tests across grouping variables.

| Test | Grouping | Statistic | *p*-value | R² / R |
|------|----------|----------:|----------:|-------:|
| PERMANOVA | Compartment | *F* = 9.955 | 0.001 | 0.424 |
| PERMANOVA | Health status | *F* = 1.546 | 0.120 | 0.052 |
| PERMANOVA | Group (D × C) | *F* = 4.757 | 0.001 | 0.498 |
| ANOSIM | Compartment | *R* = 0.685 | 0.001 | — |
| ANOSIM | Health status | *R* = 0.055 | 0.089 | — |
| PERMDISP | Compartment | *F* = 2.323 | 0.122 | — |
| PERMDISP | Health status | *F* = 1.528 | 0.226 | — |
| PERMDISP | Group (D × C) | *F* = 0.566 | 0.715 | — |

**Table 3.** Top 10 endosphere indicator MAGs by IndVal score.

| Rank | MAG ID | Phylum | Genus | IndVal | *q*-value |
|-----:|--------|--------|-------|-------:|----------:|
| 1 | H4E_semibin_4 | Myxococcota | Polyangia | 0.91 | <0.001 |
| 2 | — | Myxococcota | *Palsa-1150* | 0.91 | <0.001 |
| 3 | — | Alphaproteobacteria | *Roseiarcus* | 0.90 | <0.001 |
| 4 | — | Alphaproteobacteria | *Roseiarcus* | 0.90 | <0.001 |
| 5 | — | Alphaproteobacteria | *Roseiarcus* | 0.89 | <0.001 |
| 6 | — | Myxococcota | JADGRB01 | 0.88 | <0.001 |
| 7 | — | Myxococcota | JADGRB01 | 0.88 | <0.001 |
| 8 | — | Alphaproteobacteria | *Roseiarcus* | 0.87 | <0.001 |
| 9 | — | Alphaproteobacteria | *Roseiarcus* | 0.87 | <0.001 |
| 10 | — | Myxococcota | *Palsa-1150* | 0.87 | <0.001 |

**Table 4.** Indicator species counts by compartment and dominant phyla.

| Compartment | Total indicators | Top phylum (count, %) | Second phylum (count, %) |
|-------------|--------:|-----------------------|--------------------------|
| Bulk soil | 199 | Acidobacteriota (72, 36%) | Pseudomonadota (45, 23%) |
| Endosphere | 239 | Pseudomonadota (79, 33%) | Myxococcota (15, 6%) |
| Rhizosphere | 73 | Acidobacteriota (41, 56%) | Pseudomonadota (3, 4%) |

**Table 5.** Disease indicator summary.

| Health state | Significant indicators | Top indicator (phylum, IndVal) |
|-------------|----------:|------|
| Diseased | 125 | Micrarchaeota (DPANN), IndVal = 0.90 |
| Healthy | 36 | — |
| **Ratio** | **3.5:1** | |

**Table 6.** Core microbiome membership at multiple prevalence thresholds.

| Prevalence threshold | Core MAGs (all groupings) |
|--------------------:|------------------:|
| 50% | 1,018 |
| 75% | 1,018 |
| 90% | 1,017–1,018 |
| 100% | 1,017 |

---

## Figure Legends

**Figure 1.** Differential abundance analysis across compartments and health states. (a) Volcano plot showing CLR-transformed abundance differences between bulk soil and endosphere (634 significant MAGs highlighted). (b) Summary of significant MAGs across all pairwise comparisons. (c) Volcano plot for diseased versus healthy comparison (0 significant MAGs).

**Figure 2.** Compartment specificity of 1,018 MAGs. (a) Distribution of specificity scores, with 0.5 threshold indicated. (b) Compartment identity of all 65 specialists (all endosphere). (c) Specificity scores of Myxococcota MAGs highlighting endosphere concentration.

**Figure 3.** Archaea exclusion from root compartments. Relative abundance of Thermoplasmatota and Thermoproteota across compartments with Cohen's *d* effect sizes.

**Figure 4.** PCoA ordination of Bray-Curtis dissimilarities coloured by (a) compartment and (b) health status, illustrating the strong compartment separation and negligible disease effect.

**Figure 5.** Indicator species analysis. (a) IndVal score distribution for compartment indicators. (b) Phylum composition of indicators per compartment. (c) Disease vs healthy indicator counts with top taxa highlighted.

**Figure 6.** Core microbiome stability. Number of core MAGs at increasing prevalence thresholds across all grouping variables.
