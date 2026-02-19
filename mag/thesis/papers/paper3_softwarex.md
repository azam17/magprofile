# magprofile: A Python toolkit for genome-resolved community profiling, functional annotation, and co-occurrence network analysis of metagenome-assembled genomes

---

## Code Metadata

| Nr. | Code metadata description | Value |
|-----|--------------------------|-------|
| C1  | Current code version | 0.4.3 |
| C2  | Permanent link to code/repository used for this code version | https://github.com/[username]/magprofile |
| C3  | Permanent link to reproducible capsule | *To be provided* |
| C4  | Legal code license | MIT |
| C5  | Code versioning system used | git |
| C6  | Software code languages, tools, and services used | Python 3.10+ |
| C7  | Compilation requirements, operating environments, and dependencies | NumPy ≥ 1.24, SciPy ≥ 1.10, NetworkX, Matplotlib, Seaborn, Click ≥ 8.0, Pandas |
| C8  | If available, link to developer documentation/manual | https://github.com/[username]/magprofile/wiki |
| C9  | Support email for questions | *Corresponding author email* |

---

## Abstract

Metagenome-assembled genomes (MAGs) have become central to culture-independent microbiology, yet analysing MAG abundance tables across samples, functional annotations, and inter-taxon interactions typically requires stitching together disparate R packages and custom scripts with no unified interface. We present **magprofile**, an open-source Python toolkit that integrates three analytical pillars — community ecology, functional profiling, and co-occurrence network analysis — into a single command-line application. The ecology module (`magprofile report`) computes alpha and beta diversity, PERMANOVA with PERMDISP validation, indicator species analysis, core microbiome identification, and compartment specificity scoring. The functional module (`magprofile func-report`) parses DRAM annotations to construct KEGG Orthology, CAZyme, and plant growth-promoting rhizobacteria (PGPR) trait matrices, then performs differential pathway analysis and functional redundancy assessment. The network module (`magprofile net-report`) builds compositionally-aware co-occurrence networks using proportionality (phi), applies degree-preserving null models, identifies keystone taxa and hub-bridge roles via the Guimerà–Amaral framework, and supports differential network comparison between experimental groups. All modules produce publication-ready figures and self-contained HTML reports. The package comprises 6,090 lines of source code across 20 modules, with 158 unit tests (2,670 lines) achieving comprehensive coverage. magprofile was developed for and validated on a 30-sample, 1,018-MAG oil palm root microbiome dataset, where it revealed that root compartment — not Ganoderma basal stem rot disease status — is the dominant driver of community assembly (PERMANOVA R² = 0.424), functional differentiation (2,000+ differential pathways), and network topology. magprofile is pip-installable, requires no compiled dependencies, and is freely available under the MIT license.

**Keywords:** metagenome-assembled genomes; community profiling; co-occurrence networks; functional annotation; microbiome; Python

---

## 1. Motivation and significance

Shotgun metagenomics coupled with genome binning now routinely yields hundreds to thousands of metagenome-assembled genomes (MAGs) per study [1,2]. Once MAGs are recovered and taxonomically classified, researchers face a multi-layered analytical challenge: quantifying community composition across samples, annotating functional capacity, and inferring ecological interactions among taxa. Each of these tasks is typically addressed by separate software tools originating from different language ecosystems.

For community ecology, R packages such as vegan [3] and phyloseq [4] dominate, but they operate on operational taxonomic unit (OTU) or amplicon sequence variant (ASV) tables and require manual adaptation for MAG-level data. Functional profiling tools such as DRAM [5] generate raw annotation tables but provide no statistical framework for comparing functional repertoires across experimental groups. Network inference tools such as SparCC [6], FlashWeave [7], and SpiecEasi [8] focus on correlation estimation but lack integrated downstream analysis — keystone identification, null model validation, and hub-bridge classification must be scripted separately.

This fragmentation imposes three practical costs on researchers:

1. **Language switching.** Analyses commonly span Python (for upstream bioinformatics), R (for statistics), and Bash (for glue scripts), creating reproducibility barriers and increasing the likelihood of data transformation errors between tools.

2. **Compositional data violations.** Many existing tools apply Pearson or Spearman correlation to relative abundance data without compositional correction, generating spurious associations [9,10]. Researchers must independently implement centered log-ratio (CLR) transformation for differential abundance and proportionality metrics for network construction.

3. **Lack of integration.** Ecological findings (e.g., differential taxa), functional annotations (e.g., PGPR traits), and network roles (e.g., keystone status) exist in separate output files with no framework for cross-referencing. Identifying whether a keystone taxon carries plant-beneficial genes, for instance, requires ad hoc scripting.

magprofile addresses these gaps by providing a unified Python toolkit that handles all three analytical pillars through a consistent command-line interface, uses compositionally-aware statistical methods throughout, and produces integrated cross-analysis outputs. It specifically targets the MAG-centric workflow — accepting MAG × sample abundance matrices, MAG-level taxonomy, and DRAM functional annotations as inputs — filling a niche between upstream genome recovery tools (e.g., MetaWRAP [11], ATLAS [12]) and downstream visualisation platforms.

The software was developed in the context of an oil palm (*Elaeis guineensis*) root microbiome study investigating the impact of Ganoderma boninense basal stem rot (BSR) on soil microbial communities. The three-module design directly reflects the biological questions arising from this system: who is there (ecology), what can they do (function), and how do they interact (network). However, the toolkit is domain-agnostic and applicable to any MAG-based study with sample metadata and optional functional annotations.

## 2. Software description

### 2.1. Software architecture

magprofile is organised as a Python package of 20 modules (6,090 lines of source code) structured into three analytical pillars plus shared infrastructure (Fig. 1):

- **Core I/O** (`io.py`): Defines data types (`AbundanceTable`, `TaxonomyTable`, `SampleMetadata`) and loaders for tab-separated MAG abundance matrices, GTDB-style taxonomy, and sample metadata.
- **Ecology pillar** (`diversity.py`, `beta.py`, `ordination.py`, `differential.py`, `indicator.py`, `core_microbiome.py`, `compartment_specificity.py`, `stats.py`): Eight modules covering alpha diversity (Shannon, Simpson, richness, Pielou evenness), beta diversity (Bray–Curtis, Jaccard), ordination (PCoA, NMDS), compositional differential abundance (CLR + Welch's t-test with BH-FDR), indicator species (IndVal with permutation testing), core microbiome at multiple prevalence thresholds, information-theoretic compartment specificity, and permutational statistics (PERMANOVA, ANOSIM, PERMDISP, pairwise PERMANOVA, variance partitioning, Mantel test).
- **Functional pillar** (`func_io.py`, `func_profile.py`, `func_report.py`): Three modules that parse DRAM annotation tables into KEGG Orthology (KO), CAZyme, and PGPR marker gene matrices; project function × MAG matrices to function × sample abundance; compute differential pathway abundance between groups; assess functional redundancy as Shannon diversity of function carriers; and perform keystone–PGPR cross-referencing via Fisher's exact test.
- **Network pillar** (`net_correlation.py`, `net_topology.py`, `net_report.py`): Three modules implementing proportionality (phi) for compositionally-aware co-occurrence estimation [9], threshold sensitivity analysis across phi percentiles, graph topology metrics (degree, betweenness, closeness, modularity via Louvain), keystone taxa identification (high centrality + low abundance), degree-preserving null models for modularity validation, hub-bridge classification via the Guimerà–Amaral z-P framework [13], and differential network comparison between experimental groups.
- **Visualisation and reporting** (`plots.py`, `html_report.py`): 15 plotting functions producing publication-quality figures (boxplots, ordination scatter with 95% confidence ellipses, volcano plots, heatmaps, network layouts, hub-bridge scatter) and self-contained HTML reports with embedded figures and statistical tables.
- **CLI** (`cli.py`): Click-based command-line interface exposing seven subcommands: `diversity`, `ordination`, `differential`, `indicator`, `report` (full ecology pipeline), `func-report` (full functional pipeline), and `net-report` (full network pipeline).

All numerical computation uses NumPy and SciPy with no compiled C extensions, ensuring straightforward installation on any platform with Python ≥ 3.10.

### 2.2. Software functionalities

**Table 1.** Summary of magprofile analytical capabilities.

| Module | Command | Key outputs |
|--------|---------|-------------|
| Ecology | `magprofile report` | Alpha diversity (4 metrics), beta diversity (Bray–Curtis, Jaccard), PCoA/NMDS ordination, PERMANOVA/ANOSIM/PERMDISP, pairwise PERMANOVA with FDR, variance partitioning, differential abundance (CLR + Cohen's d), indicator species (IndVal), core microbiome, compartment specificity, HTML report |
| Functional | `magprofile func-report` | KO abundance matrix, CAZyme class summary, PGPR trait enrichment (chi-square/Fisher + FDR), differential pathway abundance, functional redundancy (Shannon of carriers), keystone–PGPR cross-reference, HTML report |
| Network | `magprofile net-report` | Phi proportionality matrix, global/group networks, threshold sensitivity, topology metrics (degree, betweenness, closeness, modularity), keystone taxa, null model z-score, module composition, hub-bridge classification (z-P), differential networks, HTML report |

**Compositional data handling.** Microbiome abundance data are inherently compositional — relative abundances sum to a constant — violating assumptions of standard parametric tests and correlation metrics [10]. magprofile addresses this at two levels: (1) differential abundance uses CLR transformation with pseudocount addition before Welch's t-test, and (2) co-occurrence network construction uses the phi proportionality metric [9], which quantifies the log-ratio variance between taxa pairs rather than correlation, avoiding the spurious associations generated by Pearson or Spearman on compositional data.

**Statistical rigour.** All hypothesis tests that support permutation (PERMANOVA, ANOSIM, IndVal, Mantel) use permutation-based p-values rather than parametric approximations. Multiple comparison correction uses Benjamini–Hochberg FDR throughout. Effect sizes (Cohen's d) are reported alongside p-values for differential abundance. Network modularity is validated against degree-preserving null models with z-score reporting.

**Keystone taxa identification.** Keystone taxa are defined as MAGs with disproportionately high network centrality (betweenness and/or closeness) relative to their abundance, following the ecological keystone concept [14]. The algorithm ranks MAGs by a composite score combining normalised betweenness centrality and inverse normalised abundance, identifying taxa that are rare but structurally critical to the co-occurrence network.

**Hub-bridge classification.** Using the Guimerà–Amaral z-P framework [13], magprofile classifies each node by its within-module degree z-score (hub: z > 2.5) and participation coefficient P (bridge: P > 0.62). This distinguishes module hubs (strongly connected within their own module) from connector nodes (linking multiple modules), providing ecological insight into community organisation.

### 2.3. Input and output

**Inputs:**
- MAG × sample abundance table (TSV, tab-delimited; rows = MAGs, columns = samples)
- Sample metadata (TSV; must include a grouping variable column, e.g., "compartment")
- Taxonomy table (optional TSV; GTDB-style with columns: domain, phylum, class, order, family, genus, species)
- DRAM annotations (TSV; required for `func-report` only)

**Outputs (per module):**
- CSV files for all numerical results (diversity metrics, differential abundance, network topology, etc.)
- PDF figures for all visualisations
- Self-contained HTML report with embedded figures and summary statistics

### 2.4. Installation

```bash
pip install magprofile
```

Or from source:

```bash
git clone https://github.com/[username]/magprofile.git
cd magprofile
pip install -e ".[dev]"
```

Dependencies are resolved automatically via pip. No system-level compilation is required.

## 3. Illustrative examples

### 3.1. Full ecology analysis

We demonstrate magprofile on a dataset of 1,018 MAGs (567 high-quality, 451 medium-quality) recovered from 30 oil palm root zone soil samples spanning three compartments (bulk soil, rhizosphere, endosphere) and two health states (healthy, Ganoderma BSR-diseased), with five biological replicates per group.

```bash
magprofile report \
    --abundance abundance_table.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --group compartment \
    --permutations 999 \
    --output ecology_results/
```

This single command executes the complete ecology pipeline: alpha diversity computation (Shannon H = 5.994 ± 0.227, Simpson D = 0.994 ± 0.002 across 30 samples), Bray–Curtis beta diversity, PCoA ordination, PERMANOVA (compartment: R² = 0.424, p = 0.001; health: R² = 0.052, p = 0.120), PERMDISP validation (all p > 0.12, confirming homogeneous dispersions), ANOSIM (compartment: R = 0.685, p = 0.001), pairwise PERMANOVA with FDR correction (bulk vs endosphere: R² = 0.462; endosphere vs rhizosphere: R² = 0.331; bulk vs rhizosphere: R² = 0.194), differential abundance (634 significant MAGs between bulk and endosphere at FDR < 0.05), indicator species analysis (199 bulk, 239 endosphere, 73 rhizosphere indicators), core microbiome (1,017 of 1,018 MAGs present at 100% prevalence), and compartment specificity (65 specialists, all endosphere-associated). All results, figures, and an HTML report are written to `ecology_results/`.

### 3.2. Functional profiling

```bash
magprofile func-report \
    --abundance abundance_table.tsv \
    --dram-annotations dram_annotations.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --group compartment \
    --output func_results/
```

This analyses 567 DRAM-annotated HQ MAGs, constructing KO (6,233 functions), CAZyme (532 families), and PGPR (13 traits) matrices. Key findings include: nitrogen fixation genes (nifH, nifK) significantly enriched in the endosphere (FDR < 0.05), ACC deaminase (acdS) enriched in the endosphere (15 vs 5 bulk MAGs, FDR < 0.01), functional redundancy significantly higher in root compartments than bulk soil (Shannon 2.99 vs 2.71, Kruskal–Wallis p < 0.001), and 2,000+ differentially abundant pathways between bulk and endosphere but zero between diseased and healthy samples.

### 3.3. Network analysis

```bash
magprofile net-report \
    --abundance abundance_table.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --group compartment \
    --threshold 5.0 \
    --threshold-mode global \
    --min-prevalence 0.5 \
    --null-iterations 1000 \
    --null-jobs 4 \
    --output net_results/
```

The network module constructs a global co-occurrence network (25,883 edges, density = 0.050) from all 1,018 MAGs, validates modularity (0.226) against 1,000 degree-preserving null models (z = 81.36, p < 0.001), identifies 258 modules with the largest dominated by Acidobacteriota (341 MAGs), and classifies keystone taxa (top: *Terracidiphilus* sp., Acidobacteriota, keystone score = 0.981). Group-level networks reveal 27% higher density in diseased samples (0.117 vs 0.092), indicating structural dysbiosis without compositional change.

### 3.4. Cross-pillar integration

The three modules produce outputs in consistent CSV formats that enable cross-referencing. For example, keystone taxa identified by the network module can be matched against PGPR carriers from the functional module:

```python
import pandas as pd

keystones = pd.read_csv("net_results/keystones_global.csv")
pgpr = pd.read_csv("func_results/pgpr_enrichment.csv")

# Identify keystones carrying plant-beneficial traits
keystone_ids = set(keystones["mag_id"])
pgpr_carriers = set(pgpr.query("count > 0")["mag_id"])
beneficial_keystones = keystone_ids & pgpr_carriers
```

In the oil palm dataset, this cross-referencing revealed that Acidobacteriota keystones carry phosphate-solubilising genes (gcd), connecting network structural importance to functional plant-growth promotion.

## 4. Impact

### 4.1. Current use

magprofile was developed for and applied to the analysis of a 30-sample oil palm root microbiome dataset as part of doctoral research at Universiti Putra Malaysia. The three-pillar analysis (ecology, function, network) generated results for two peer-reviewed manuscripts: (1) a community ecology paper revealing compartment-driven assembly under Ganoderma BSR, and (2) a functional and network paper demonstrating the three-layer insurance hypothesis — taxonomic redundancy, functional backup, and modular network architecture collectively buffer the microbiome against disease perturbation.

### 4.2. Comparison with existing tools

**Table 2.** Feature comparison of magprofile with existing microbiome analysis tools.

| Feature | magprofile | vegan (R) | phyloseq (R) | QIIME 2 | FlashWeave | SparCC |
|---------|-----------|-----------|--------------|---------|------------|--------|
| MAG-native input | Yes | No | No | No | No | No |
| Alpha/beta diversity | Yes | Yes | Yes | Yes | No | No |
| PERMANOVA + PERMDISP | Yes | Yes | No | Yes | No | No |
| Compositional differential | CLR + Welch | Manual | DESeq2 wrap | ANCOM-BC | No | No |
| Indicator species (IndVal) | Yes | Yes | No | No | No | No |
| DRAM functional parsing | Yes | No | No | No | No | No |
| PGPR trait enrichment | Yes | No | No | No | No | No |
| Functional redundancy | Yes | No | No | No | No | No |
| Compositional network (phi) | Yes | No | No | No | Partial | Yes |
| Null model validation | Yes | No | No | No | No | No |
| Keystone identification | Yes | No | No | No | No | No |
| Hub-bridge (z-P) | Yes | No | No | No | No | No |
| Integrated HTML report | Yes | No | No | Yes | No | No |
| Single-language | Python | R | R | Python/R | Julia | Python |
| CLI interface | Yes | No | No | Yes | Yes | No |

magprofile is, to our knowledge, the first tool to combine compositional community ecology, DRAM-based functional profiling, and compositionally-aware co-occurrence network analysis with null model validation in a single Python package.

### 4.3. Software quality

The package includes 158 unit tests (2,670 lines of test code) covering all 20 modules. Tests validate numerical correctness (e.g., Shannon diversity against hand-calculated values), edge cases (e.g., zero-abundance samples, single-MAG communities), statistical properties (e.g., PERMANOVA p-value distribution under null), and integration (e.g., full pipeline execution with synthetic data). The test suite executes in under 10 seconds on a standard laptop.

### 4.4. Potential impact

magprofile enables several research directions that are currently difficult with fragmented tools:

1. **Integrated multi-pillar microbiome analysis.** By unifying ecology, function, and network analysis with consistent data structures, researchers can ask cross-cutting questions (e.g., "Are keystone taxa functionally distinct?") without ad hoc data merging.

2. **Compositionally-correct workflows by default.** CLR transformation for differential abundance and phi proportionality for networks are applied automatically, reducing the risk of spurious findings from inappropriate statistical methods.

3. **MAG-centric workflows.** Unlike amplicon-oriented tools, magprofile operates natively on MAG abundance tables and integrates DRAM annotations, reflecting the growing adoption of genome-resolved metagenomics.

4. **Reproducible reporting.** Self-contained HTML reports embed all figures, tables, and statistical summaries, facilitating peer review and supplementary material generation.

## 5. Conclusions

magprofile provides a unified, compositionally-aware Python toolkit for the three core analytical tasks in MAG-based microbiome research: community ecology, functional profiling, and co-occurrence network analysis. Its design fills a practical gap between upstream genome recovery tools and downstream statistical environments by providing domain-specific analyses (IndVal, phi proportionality, null model validation, PGPR enrichment, hub-bridge classification) through a simple command-line interface. The package is freely available under the MIT license and is installable via pip with no compiled dependencies.

---

## Acknowledgements

This work was supported by Universiti Putra Malaysia. We thank the Malaysian Palm Oil Board (MPOB) for providing access to oil palm plantation sampling sites.

---

## References

[1] Parks DH, Rinke C, Chuvochina M, Chaumeil PA, Woodcroft BJ, Evans PN, Hugenholtz P, Tyson GW. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. Nature Microbiology. 2017;2(11):1533–1542.

[2] Nayfach S, Roux S, Seshadri R, Udwary D, Varghese N, Schulz F, Wu D, Paez-Espino D, Chen IM, Huntemann M, Palaniappan K. A genomic catalog of Earth's microbiomes. Nature Biotechnology. 2021;39(4):499–509.

[3] Oksanen J, Simpson GL, Blanchet FG, Kindt R, Legendre P, Minchin PR, et al. vegan: Community Ecology Package. R package version 2.6-4. 2022.

[4] McMurdie PJ, Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;8(4):e61217.

[5] Shaffer M, Borton MA, McGivern BB, Zayed AA, La Rosa SL, Solden LM, Liu P, Narrowe AB, Rodríguez-Ramos J, Bolduc B, Gazitúa MC. DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic Acids Research. 2020;48(16):8883–8900.

[6] Friedman J, Alm EJ. Inferring correlation networks from genomic survey data. PLoS Computational Biology. 2012;8(9):e1002687.

[7] Tackmann J, Rodrigues JFM, von Mering C. Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell Systems. 2019;9(3):286–296.

[8] Kurtz ZD, Müller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA. Sparse and compositionally robust inference of microbial ecological networks. PLoS Computational Biology. 2015;11(5):e1004226.

[9] Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S, Bähler J. Proportionality: a valid alternative to correlation for relative data. PLoS Computational Biology. 2015;11(3):e1004075.

[10] Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology. 2017;8:2224.

[11] Uritskiy GV, DiRuggiero J, Taylor J. MetaWRAP — a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome. 2018;6(1):158.

[12] Kieser S, Brown J, Zdobnov EM, Berthel M, Mangeat B. ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data. BMC Bioinformatics. 2020;21(1):257.

[13] Guimerà R, Amaral LAN. Functional cartography of metabolic networks. Nature. 2005;433(7028):895–900.

[14] Berry D, Widder S. Deciphering microbial interactions and detecting keystone species with co-occurrence networks. Frontiers in Microbiology. 2014;5:219.

[15] Anderson MJ. A new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32–46.

[16] Anderson MJ. Distance-based tests for homogeneity of multivariate dispersions. Biometrics. 2006;62(1):245–253.

[17] Clarke KR. Non-parametric multivariate analyses of changes in community structure. Australian Journal of Ecology. 1993;18(1):117–143.

[18] Dufrêne M, Legendre P. Species assemblages and indicator species: the need for a flexible asymmetrical approach. Ecological Monographs. 1997;67(3):345–366.

[19] Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs. 1957;27(4):325–349.

[20] Allison SD, Martiny JBH. Resistance, resilience, and redundancy in microbial communities. Proceedings of the National Academy of Sciences. 2008;105(Suppl 1):11512–11519.

---

## Figure legends

**Fig. 1.** Software architecture of magprofile. The package is organised into three analytical pillars (ecology, function, network) sharing a common I/O layer and producing integrated outputs via a CLI interface. Arrows indicate data flow between modules.

**Fig. 2.** Illustrative output from the ecology module. (A) PCoA ordination of 30 oil palm root zone samples coloured by compartment, showing clear separation of bulk soil from root-associated communities. (B) Box plots of Shannon diversity by compartment. (C) Volcano plot of differential abundance (bulk vs endosphere), with 634 significant MAGs highlighted.

**Fig. 3.** Illustrative output from the network module. (A) Global co-occurrence network (1,018 nodes, 25,883 edges) coloured by phylum, showing the dominant Acidobacteriota module. (B) Hub-bridge scatter plot (within-module degree z vs participation coefficient P) classifying MAGs into network roles. (C) Threshold sensitivity analysis showing modularity stability across phi percentile thresholds.

**Fig. 4.** Cross-pillar integration example. Keystone taxa identified by the network module (red nodes) are cross-referenced with PGPR carriers from the functional module (green ring), revealing that structurally important taxa carry plant-beneficial genes.
