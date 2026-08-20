# FRASER-visualizations
This repository contains the code to create different benchmarking visualizations on pyfraser

## Usage

### `Proportion of splicing outliers with rare splice-disrupting variants`

Compares an arbitrary set of FRASER runs (any mix of Python/R outputs) on the same cohort and plots splicing-outlier enrichment for supporting variants.

```python
from Proportion_with_variants import compare_splicing_runs

runs = [
    {"path": "results/fraser_summary_all_junctions.csv.gz", "kind": "python", "label": "Python"}, # the csv output containing all junctions
    {"path": "results/raw-local-Spleen_results_all.tsv.gz", "kind": "r", "label": "R"}, # the tsv output containing all junctions
]

compare_splicing_runs(
    runs=runs,
    vep_path="data/vep_annotated_variants.csv",
    out_path="plots/splicing_enrichment.png",
    ylabel="Proportion of splicing outliers with\nrare splice-disrupting variants",
    outlier_source="Python",  # must match one run's "label". The number of outliers from this source will draw a vertical dashed line
)
```

### `Venn diagram of overlapping gene/junction sets`

Compares the gene/junction sets of two FRASER outputs (any mix of Python/R) and draws a Venn diagram of the overlap.

```python
from Venn_diagram import compare_outputs

compare_outputs(
    source_1="results/fraser_summary_filtered_junctions.csv",  # path to a CSV/TSV, or an already-loaded DataFrame containing only the filtered outliers
    source_2="results/raw-local-Spleen_results_per_junction.tsv", # similar to source_1
    level="junction",  # "junction" or "gene" level comparison
    cohort="Spleen", # the name of the cohort for the title
    name_1="Python", # label of the first source
    name_2="R", # label of the second source
    out_path="plots/venn_junction_overlap.png",  # if given, the figure is saved here
)
```
