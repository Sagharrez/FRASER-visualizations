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

### `Step timing comparison across FRASER runs`

Compares an arbitrary set of FRASER runs (any mix of Python/R outputs) on aligned pipeline-step buckets and plots per-cohort timing bars. It plots timings for all the cohorts as subfolders.

```python
from Compare_timing import compare_timing, plot_timing

runs = [
    Run("R",          "R",      "R_V2/AE_CPU_v2/{cohort}/*_timing.csv"),
    Run("Python CPU", "Python", "Python/{cohort}/AE_16ep_CPU/protrider.log"),
    Run("Python GPU", "Python", "Python/{cohort}/AE_16ep/protrider.log"),
]
fig = plot_timing(
    compare_timing(runs), 
    runs, 
    baseline="R", 
    title="FRASER step timing by cohort on autoencode with 16 epochs", 
    output_path = 'plots/fraser_timing_R_Py.png')
```

### `Subset recovery across tissues`

For each tissue's own outlier (sample, junction) pairs, checks how many are also flagged as outliers in a joint (Subset/union) FRASER run, and plots a horizontal stacked bar of recovered vs. missed pairs per tissue.

```python
from Subset_recovery import plot_subset_recovery

root = "/path/to/root"
fig, ax, summary = plot_subset_recovery(
    f"{root}/Subset/PCA/fraser_summary_filtered_junctions.csv", # path to early fusion run filtered junctions
    {
        "Brain": f"{root}/Brain/PCA/fraser_summary_filtered_junctions.csv", # paths to each of the cohorts' run filtered junctions
        "Lung": f"{root}/Lung/PCA/fraser_summary_filtered_junctions.csv",
        "Pancreas": f"{root}/Pancreas/PCA/fraser_summary_filtered_junctions.csv",
        "Spleen": f"{root}/Spleen/PCA/fraser_summary_filtered_junctions.csv",
    },
)
fig.savefig("plots/subset_recovery.png", dpi=150, bbox_inches="tight")
```
