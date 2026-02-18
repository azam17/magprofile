# magfunc + magnet Design Document

**Date:** 2026-02-18
**Status:** Approved
**Scope:** Add functional profiling (Chapter 2) and network analysis (Chapter 3) to the existing mag/ package

---

## Context

magprofile v0.2.1 (Chapter 1: "Who's there and where?") revealed:
- Myxococcota predatory guild in endosphere (15 MAGs, 10x enrichment)
- Roseiarcus photoheterotrophs inside opaque roots
- Acidobacteriota rhizosphere dominance (56% of indicators)
- Near-complete archaea exclusion from roots (Cohen's d > 6.0)
- 3 disease-associated Verrucomicrobiota MAGs
- Ultrasmall organisms (CPR/DPANN) as dysbiosis markers

These findings demand functional annotation ("what genes enable this?") and interaction analysis ("who are they interacting with?").

## Architecture Decision

Both tools are added as modules within the existing `mag/` package, sharing data types (AbundanceTable, TaxonomyTable, SampleMetadata) and utilities (CLR transform, BH-FDR, permutation tests).

## Input Data

- **magfunc:** DRAM output (`annotations.tsv` + `distill/` directory) + existing AbundanceTable + Metadata
- **magnet:** Existing AbundanceTable + Metadata + (optionally) FunctionalTable from magfunc

---

## Unified Data Model

### FunctionalTable (new, in func_io.py)

```python
@dataclass(frozen=True)
class DRAMAnnotation:
    gene_id: str
    mag_id: str
    ko_id: str | None
    kegg_hit: str | None
    cazy_id: str | None
    pfam_id: str | None

@dataclass(frozen=True)
class FunctionalTable:
    function_ids: list[str]     # KO IDs, module IDs, or CAZy families
    mag_ids: list[str]
    values: np.ndarray          # n_functions x n_mags (binary, counts, or completeness 0-1)
    function_type: str          # "ko" | "module" | "cazy" | "pgpr"

    def to_sample_abundance(self, abundance: AbundanceTable) -> AbundanceTable:
        """Multiply function x MAG @ MAG x sample -> function x sample abundance."""
```

### Network Types (new, in net_correlation.py / net_topology.py)

```python
@dataclass(frozen=True)
class NetworkResult:
    mag_ids: list[str]
    adjacency: np.ndarray       # n_mags x n_mags correlation values
    edges: list[tuple[str, str, float]]  # (mag1, mag2, weight) filtered by threshold
    threshold: float

@dataclass(frozen=True)
class NetworkTopology:
    mag_ids: list[str]
    degree: np.ndarray
    betweenness: np.ndarray
    closeness: np.ndarray
    modularity: float
    module_assignments: np.ndarray
    hub_scores: np.ndarray

@dataclass(frozen=True)
class KeystoneTaxaResult:
    mag_ids: list[str]
    keystone_scores: np.ndarray
    is_keystone: np.ndarray     # boolean
    metrics: dict[str, np.ndarray]
```

---

## magfunc — Functional Profiling (Chapter 2)

### Module: func_io.py — DRAM Parsing

- `load_dram_annotations(path) -> list[DRAMAnnotation]` — parse annotations.tsv
- `load_dram_distill(distill_dir) -> FunctionalTable` — parse module completeness from distill/
- `build_functional_table(annotations, function_type) -> FunctionalTable` — aggregate annotations into function x MAG matrix
- Built-in PGPR marker dict: nifH (K02588), nifD (K02586), nifK (K02591), pqqC (K06137), gcd (K00111), acdS (K01505), ipdC (K04103), entA-F, etc.

### Module: func_profile.py — Analysis

- `pathway_abundance(func_table, abundance) -> AbundanceTable` — function x sample matrix weighted by MAG abundance
- `differential_pathway(func_table, abundance, metadata, grouping_var, g1, g2) -> DifferentialAbundanceResult` — reuses existing differential_abundance() pipeline
- `functional_redundancy(func_table, abundance, metadata, grouping_var) -> dict` — per function per group: number of carrier MAGs + Shannon diversity of carriers
- `pgpr_profile(annotations) -> FunctionalTable` — subset to curated PGPR markers, binary per MAG
- `cazyme_profile(annotations) -> FunctionalTable` — CAZy families (GH/GT/PL/CE/AA), counts per MAG

### Module: func_report.py — Output

- `generate_func_report(abundance, func_table, taxonomy, metadata, grouping_var, output_dir)`
- CSVs: pathway_abundance.csv, pathway_differential_*.csv, functional_redundancy.csv, pgpr_traits.csv, cazyme_summary.csv
- Plots: pathway heatmap (top differential pathways), PGPR binary heatmap, CAZyme class barplot per group, redundancy barplot
- HTML: func_report.html (same self-contained style as magprofile)

---

## magnet — Network Analysis (Chapter 3)

### Module: net_correlation.py — Compositional Correlation

- `proportionality(abundance, metric="phi") -> np.ndarray` — phi = var(log(x/y)), compositionally aware
- `build_network(corr_matrix, threshold, min_prevalence) -> NetworkResult` — threshold edges, filter rare MAGs
- Default: keep MAGs in >= 50% of samples, threshold phi at 5th percentile (strong associations)

### Module: net_topology.py — Network Metrics

- `compute_topology(network) -> NetworkTopology` — degree, betweenness, closeness, Louvain modularity (via networkx)
- `identify_keystones(topology, abundance) -> KeystoneTaxaResult` — high betweenness + high degree + low abundance = keystone
- `differential_network(net1, net2) -> dict` — edges gained/lost between two compartment networks

### Module: net_report.py — Output

- `generate_net_report(abundance, taxonomy, metadata, grouping_var, output_dir)`
- CSVs: network_edges.csv, network_topology.csv, keystone_taxa.csv, differential_network_*.csv
- Plots: network graph (force-directed, nodes by phylum), degree distribution, modularity barplot
- HTML: net_report.html

---

## CLI Integration (cli.py additions)

```
magprofile func-report  --abundance --dram-annotations --dram-distill --taxonomy --metadata --group --output
magprofile net-report   --abundance --taxonomy --metadata --group --threshold --min-prevalence --output
```

---

## Dependencies

- magfunc: No new dependencies (DRAM output is TSV, uses existing numpy/scipy)
- magnet: `networkx` (graph algorithms: Louvain modularity, betweenness centrality, force-directed layout)

---

## File Layout

```
mag/
├── func_io.py              # DRAM parsing -> FunctionalTable
├── func_profile.py         # pathway abundance, differential, redundancy, PGPR, CAZyme
├── func_report.py          # CSV + plot + HTML for functional results
├── net_correlation.py      # proportionality correlation, network construction
├── net_topology.py         # modularity, centrality, keystones, differential networks
├── net_report.py           # CSV + plot + HTML for network results
├── cli.py                  # ADD func-report + net-report commands
├── plots.py                # ADD pathway_heatmap, pgpr_matrix, cazyme_bar, network_graph
└── tests/
    ├── test_func_io.py     # DRAM parsing, FunctionalTable construction
    ├── test_func_profile.py # pathway abundance, redundancy, PGPR, CAZyme
    ├── test_func_report.py # end-to-end functional report
    ├── test_net_correlation.py # proportionality, network construction
    ├── test_net_topology.py    # topology metrics, keystones
    └── test_net_report.py      # end-to-end network report
```

---

## Testing Strategy

- Synthetic fixtures: extend existing fixtures.py with functional annotations (random KO/CAZy assignments to existing 20 MAGs)
- Unit tests per module (same pattern as existing tests)
- Integration test: full report generation with synthetic data
- Verification: run on real MPOB data after implementation

## Out of Scope

- Beta-NTI / RC_bray assembly analysis (requires phylogenetic tree)
- HGT detection (requires whole-genome alignment)
- Running DRAM annotation (user provides pre-computed output)
