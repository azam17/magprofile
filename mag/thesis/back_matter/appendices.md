# Back Matter

---

## Appendix A: Supplementary Tables

### Table A1. Complete list of endosphere indicator MAGs (n = 239)

[To be generated from magprofile eco-report output. Columns: MAG ID, Phylum, Class, Order, Family, Genus, IndVal score, p-value, q-value, Specificity score, Mean abundance per compartment]

### Table A2. Complete list of rhizosphere indicator MAGs (n = 73)

[Same format as Table A1]

### Table A3. Complete list of bulk soil indicator MAGs (n = 199)

[Same format as Table A1]

### Table A4. Complete list of disease indicator MAGs (n = 125)

[Same format as Table A1, with columns for diseased and healthy mean abundance]

### Table A5. Complete list of healthy indicator MAGs (n = 36)

[Same format as Table A4]

### Table A6. All pairwise PERMANOVA comparisons (15 group pairs)

[Columns: Comparison, Pseudo-F, p-value, q-value, R², Number of differential MAGs]

### Table A7. All pairwise differential abundance results — phylum-level summary

[Columns: Comparison, Phylum, Number of significant MAGs enriched, Number depleted, Mean effect size (Cohen's d)]

### Table A8. Complete PGPR trait carrier list

[Columns: Trait, Gene, MAG ID, Phylum, Class, Genus, Compartment with highest abundance]

### Table A9. Complete keystone taxa list (top 50)

[Columns: Rank, MAG ID, Phylum, Class, Genus, Keystone score, Betweenness centrality, Closeness centrality, Degree centrality, Mean abundance, Module membership]

### Table A10. Module composition — all 258 modules

[Columns: Module ID, Number of MAGs, Top 3 phyla with fractions, Hub count, Bridge count]

### Table A11. Differential pathway list — top 100 bulk vs endosphere

[Columns: KO ID, Pathway name, CLR mean bulk, CLR mean endo, Mann-Whitney U, p-value, q-value, Direction]

### Table A12. Differential CAZyme list — all significant families

[Columns: CAZyme family, Class, Mean abundance bulk/endo/rhizo, p-value, q-value, Direction]

---

## Appendix B: Supplementary Figures

### Figure B1. MAG quality distribution
Histogram of completeness vs contamination for all 1,018 MAGs, with MIMAG quality thresholds indicated.

### Figure B2. Assembly statistics
(a) Contig N50 distribution across assemblies. (b) Total assembly length per sample. (c) Number of contigs per assembly.

### Figure B3. PCoA ordination — all grouping variables
Principal coordinates analysis of Bray-Curtis dissimilarity: (a) coloured by compartment, (b) coloured by health status, (c) coloured by group.

### Figure B4. PERMDISP boxplots
Within-group distance to centroid for (a) compartment, (b) health status, (c) group comparisons.

### Figure B5. Alpha diversity by all grouping variables
Boxplots of Shannon, Simpson, richness, and evenness by (a) compartment, (b) health status, (c) group.

### Figure B6. Volcano plots — all pairwise compartment comparisons
(a) Bulk vs Endosphere. (b) Bulk vs Rhizosphere. (c) Endosphere vs Rhizosphere. (d) Diseased vs Healthy.

### Figure B7. Phylum composition barplots
Relative abundance of top 10 phyla across (a) compartments, (b) health states, (c) all 6 groups.

### Figure B8. Compartment specificity score distribution
Histogram with specialist (>0.5) and generalist (≤0.5) thresholds.

### Figure B9. CAZyme class abundance by compartment
Stacked barplots of GH, GT, CE, AA, CBM, PL abundances across compartments.

### Figure B10. Functional redundancy distributions
Histograms of per-function Shannon diversity by (a) compartment, (b) health status.

### Figure B11. Network visualisations at multiple thresholds
Network graphs at 1st, 3rd, 5th, 10th, 20th percentile thresholds showing modularity structure.

### Figure B12. Null model distribution
Histogram of 1,000 null model modularity values with observed value (Q = 0.226) indicated.

### Figure B13. Full z-P classification plot
All 1,018 MAGs plotted in Guimerà-Amaral z-P space, coloured by phylum.

### Figure B14. Compartment-specific and health-specific network properties
(a) Network density by compartment. (b) Network density by health state. (c) Group-level network densities.

---

## Appendix C: magprofile Software Documentation

### C.1 Installation

```bash
# Clone repository
git clone https://github.com/[username]/magprofile.git
cd magprofile

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install numpy scipy matplotlib seaborn click networkx

# Install package
pip install -e .
```

### C.2 Package Structure

```
mag/
├── __init__.py          # Package initialisation
├── cli.py               # Command-line interface (Click)
├── io.py                # Data I/O (AbundanceTable, TaxonomyTable, SampleMetadata)
├── diversity.py          # Alpha diversity metrics
├── beta.py              # Beta diversity utilities
├── ordination.py         # PCoA/NMDS ordination
├── stats.py             # Statistical tests (PERMANOVA, ANOSIM, PERMDISP)
├── differential.py       # Differential abundance analysis
├── indicator.py          # Indicator species (IndVal)
├── core_microbiome.py    # Core microbiome analysis
├── compartment_specificity.py  # Compartment specificity scoring
├── report.py            # Ecological report generation
├── plots.py             # Visualisation functions
├── html_report.py       # HTML report rendering
├── func_io.py           # DRAM annotation parsing
├── func_profile.py      # Functional profiling analyses
├── func_report.py       # Functional report generation
├── net_correlation.py    # Proportionality correlation
├── net_topology.py      # Network topology and keystone analysis
├── net_report.py        # Network report generation
└── tests/               # 158 unit tests
    ├── fixtures.py
    ├── test_beta.py
    ├── test_compartment_specificity.py
    ├── test_core_microbiome.py
    ├── test_differential.py
    ├── test_diversity.py
    ├── test_func_io.py
    ├── test_func_profile.py
    ├── test_func_report.py
    ├── test_indicator.py
    ├── test_io.py
    ├── test_net_correlation.py
    ├── test_net_report.py
    ├── test_net_topology.py
    ├── test_ordination.py
    ├── test_plots.py
    ├── test_report.py
    └── test_stats.py
```

### C.3 CLI Usage

```bash
# Ecological analysis report
magprofile eco-report \
    --abundance abundance.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --grouping compartment \
    --grouping health_status \
    --grouping group \
    --output-dir eco_report/

# Functional profiling report
magprofile func-report \
    --abundance abundance.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --dram-dir dram_annotations/ \
    --grouping compartment \
    --advanced \
    --output-dir func_report/

# Network analysis report
magprofile net-report \
    --abundance abundance.tsv \
    --taxonomy taxonomy.tsv \
    --metadata metadata.tsv \
    --grouping compartment \
    --advanced \
    --output-dir net_report/
```

### C.4 Input File Formats

**abundance.tsv:** Tab-separated, samples as columns, MAGs as rows, relative abundances as values.

**taxonomy.tsv:** Tab-separated, columns: MAG_ID, Domain, Phylum, Class, Order, Family, Genus, Species.

**metadata.tsv:** Tab-separated, columns: Sample_ID, compartment, health_status, group, [additional metadata].

**dram_annotations/:** Directory containing DRAM output files (annotations.tsv, distill.tsv).

### C.5 Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prevalence-threshold` | 0.1 | Minimum fraction of samples for MAG retention |
| `--fdr-threshold` | 0.05 | FDR cutoff for significance |
| `--permutations` | 999 | Number of permutations for PERMANOVA/ANOSIM |
| `--indval-permutations` | 9999 | Permutations for indicator species |
| `--phi-percentile` | 5 | Percentile threshold for network construction |
| `--null-model-iterations` | 1000 | Number of null models for modularity test |
| `--advanced` | False | Enable advanced analyses (redundancy, keystones, etc.) |

### C.6 Test Suite

```bash
# Run all tests
python -m pytest mag/tests/ -v

# Run specific module tests
python -m pytest mag/tests/test_func_profile.py -v
python -m pytest mag/tests/test_net_topology.py -v

# Coverage report
python -m pytest mag/tests/ --cov=mag --cov-report=html
```

---

## Appendix D: Published Paper Reprints

### D.1 Paper 1

[Author names]. Compartment-driven assembly of oil palm (*Elaeis guineensis*) root microbiomes under Ganoderma basal stem rot: a metagenome-assembled genome analysis. *FEMS Microbiology Ecology*, [volume], [pages], [year].

[Insert published PDF reprint]

### D.2 Paper 2

[Author names]. Functional resilience and interaction network architecture of oil palm root-associated microbiomes under Ganoderma basal stem rot. *Pertanika Journal of Tropical Agricultural Science*, [volume], [pages], [year].

[Insert published PDF reprint]
