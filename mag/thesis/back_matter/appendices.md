# Appendices

---

## Appendix A: Supplementary Tables

### Table S1. Complete list of 125 disease indicator MAGs
[To be generated from magprofile output CSV files]
- Columns: MAG ID, Phylum, Class, Order, Family, Genus, IndVal score, p-value, q-value, Mean abundance (diseased), Mean abundance (healthy), Cohen's d

### Table S2. Complete list of 36 healthy indicator MAGs
[Same format as Table S1]

### Table S3. All pairwise PERMANOVA results (15 group comparisons)
[To be extracted from eco-report output]
- Columns: Comparison, F-statistic, p-value, q-value, R²

### Table S4. Full differential abundance results (bulk vs endosphere, 634 MAGs)
[To be generated from differential analysis output]
- Columns: MAG ID, Phylum, Genus, CLR mean (bulk), CLR mean (endo), Mann-Whitney U, p-value, q-value, Cohen's d

### Table S5. Full DRAM annotation summary per MAG (567 HQ MAGs)
[To be extracted from DRAM output]
- Columns: MAG ID, Total KOs, Total CAZymes, GH count, GT count, CE count, AA count, CBM count, PL count, PGPR traits present

### Table S6. Complete keystone taxa list (top 50)
[To be extracted from magnet keystone output]
- Columns: MAG ID, Phylum, Class, Genus, Keystone score, Betweenness, Closeness, Degree, Mean abundance, Module

### Table S7. Module composition for all 258 modules
[To be extracted from magnet module output]
- Columns: Module ID, Size (MAGs), Top 3 phyla with fractions

### Table S8. Compartment-specific PGPR carrier MAG identities
[To be extracted from magfunc output]
- For each significant PGPR trait: MAG ID, Phylum, Genus, Compartment, Mean abundance

---

## Appendix B: Supplementary Figures

### Figure S1. MAG quality distribution
- CheckM2 completeness vs contamination scatter plot for all 1,018 MAGs
- Coloured by quality tier (HQ/MQ)

### Figure S2. Rarefaction curves
- Species accumulation curves per sample
- Confirming sufficient sequencing depth

### Figure S3. Additional ordination plots
- NMDS (non-metric multidimensional scaling) complementing PCoA
- Coloured by all grouping variables

### Figure S4. Complete volcano plot panel
- All 6 compartment pairwise comparisons + disease comparison
- 7-panel figure

### Figure S5. Phylum-level abundance barplots
- Stacked barplot per sample, ordered by compartment then health
- Top 15 phyla + "Other"

### Figure S6. Alpha diversity boxplots
- Shannon, Simpson, Evenness by compartment, health, and group
- 12-panel figure

### Figure S7. Full network visualisation
- Global network with 1,018 nodes coloured by module membership
- Spring layout or force-directed

### Figure S8. Threshold sensitivity network panels
- Network visualisation at 1st, 5th, 10th, 20th percentile thresholds
- Side-by-side comparison

### Figure S9. Keystone score vs abundance scatter
- All 1,018 MAGs; keystone MAGs highlighted
- Demonstrates "rare but important" pattern

### Figure S10. Group-level network density comparison
- Barplot of density for all 6 groups
- Highlighting diseased > healthy pattern across compartments

---

## Appendix C: magprofile Software Documentation

### C.1 Installation

```bash
# Clone repository
git clone [repository URL]
cd magprofile

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -m pytest tests/ -v
```

### C.2 CLI Reference

```
Usage: magprofile [OPTIONS] COMMAND [ARGS]...

Commands:
  eco-report    Run ecological analysis and generate report
  func-report   Run functional profiling analysis and generate report
  net-report    Run network analysis and generate report
```

#### eco-report
```
magprofile eco-report \
  --abundance-table abundance.tsv \
  --taxonomy-table taxonomy.tsv \
  --metadata metadata.tsv \
  --grouping compartment \
  --output-dir eco_report/ \
  --advanced
```

#### func-report
```
magprofile func-report \
  --abundance-table abundance.tsv \
  --taxonomy-table taxonomy.tsv \
  --metadata metadata.tsv \
  --dram-dir dram_annotations/ \
  --grouping compartment \
  --output-dir func_report/ \
  --advanced
```

#### net-report
```
magprofile net-report \
  --abundance-table abundance.tsv \
  --taxonomy-table taxonomy.tsv \
  --metadata metadata.tsv \
  --grouping compartment \
  --output-dir net_report/ \
  --advanced \
  --null-model-permutations 1000
```

### C.3 Input File Formats

#### abundance.tsv
Tab-separated, MAGs as rows, samples as columns. Values: relative abundance (0–1).

#### taxonomy.tsv
Tab-separated. Columns: MAG_ID, Domain, Phylum, Class, Order, Family, Genus, Species.

#### metadata.tsv
Tab-separated. Columns: Sample_ID, compartment, health_status, group.

#### dram_annotations/
Directory containing DRAM output files:
- `annotations.tsv` — per-gene annotations
- `product.tsv` — per-MAG functional summary

### C.4 Module Architecture

| Module | File | Functions | Purpose |
|--------|------|-----------|---------|
| Core I/O | `io.py` | 8 | Data loading and validation |
| Diversity | `diversity.py` | 6 | Alpha diversity metrics |
| Beta | `beta.py` | 2 | Beta diversity utilities |
| Ordination | `ordination.py` | 3 | PCoA, NMDS |
| Statistics | `stats.py` | 12 | PERMANOVA, ANOSIM, PERMDISP |
| Differential | `differential.py` | 4 | CLR + Mann-Whitney |
| Indicator | `indicator.py` | 3 | IndVal analysis |
| Core micro. | `core_microbiome.py` | 3 | Prevalence-based core |
| Specificity | `compartment_specificity.py` | 2 | Specificity scoring |
| Plots | `plots.py` | 15+ | Visualisation |
| Report | `report.py` | 10+ | Report generation |
| HTML report | `html_report.py` | 8+ | HTML output |
| Func I/O | `func_io.py` | 6 | DRAM parsing |
| Func profile | `func_profile.py` | 10 | Functional analysis |
| Func report | `func_report.py` | 5 | Functional report |
| Net corr. | `net_correlation.py` | 4 | Phi proportionality |
| Net topology | `net_topology.py` | 12 | Network metrics |
| Net report | `net_report.py` | 5 | Network report |
| CLI | `cli.py` | 3 | Command-line interface |

### C.5 Test Suite

158 unit tests across 18 test files. Run with:

```bash
python -m pytest tests/ -v --tb=short
```

Test coverage includes: data loading, all statistical tests, analysis functions, report generation, edge cases.

### C.6 Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥1.24 | Numerical computation |
| scipy | ≥1.11 | Statistical tests |
| matplotlib | ≥3.7 | Plotting |
| seaborn | ≥0.12 | Statistical visualisation |
| click | ≥8.1 | CLI framework |
| networkx | ≥3.1 | Graph algorithms |
| scikit-learn | ≥1.3 | Ordination |

---

## Appendix D: Published Paper Reprints

### D.1 Paper 1 (Ecology)
[Author names] (2026). Compartment-driven assembly of oil palm (*Elaeis guineensis*) root microbiomes under Ganoderma basal stem rot: a metagenome-assembled genome analysis. *FEMS Microbiology Ecology*, [volume], [pages].
[Insert reprint or proof]

### D.2 Paper 2 (Function + Network)
[Author names] (2026). Functional resilience and interaction network architecture of oil palm root-associated microbiomes under Ganoderma basal stem rot. *Pertanika Journal of Tropical Agricultural Science*, [volume], [pages].
[Insert reprint or proof]
