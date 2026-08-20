# FRASER-visualizations
This repository contains the code to create different benchmarking visualizations on pyfraser

## Usage

### `compare_splicing_runs`

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
