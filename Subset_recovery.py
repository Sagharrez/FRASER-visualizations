import matplotlib.pyplot as plt
import pandas as pd

JUNCTION_COLS = ["seqnames", "start", "end", "strand"]

COLOR_RECOVERED = "#2a78d6"  # categorical slot 1 (blue)
COLOR_MISSED = "#c9c8c0"     # neutral gray
COLOR_TEXT = "#52514e"       # secondary ink
COLOR_GRID = "#e5e4df"


def _pairs(df):
    return set(map(tuple, df[["sampleID"] + JUNCTION_COLS].to_numpy()))


def plot_subset_recovery(subset_path, tissue_paths: dict):
    """Horizontal stacked bar: for each tissue's own outlier (sample, junction) pairs,
    how many are also outliers in the Subset (union) run vs missed.

    Args:
        subset_path: path to Subset's fraser_summary_filtered_junctions.csv
        tissue_paths: dict of {tissue_name: path_to_fraser_summary_filtered_junctions.csv}
        ax: optional matplotlib Axes to draw into

    Returns:
        (fig, ax, summary_df) — summary_df has one row per tissue with
        n_total / n_recovered / n_missed / frac_recovered
    """
    subset_pairs = _pairs(pd.read_csv(subset_path))

    rows = []
    for tissue, path in tissue_paths.items():
        tissue_pairs = _pairs(pd.read_csv(path))
        n_recovered = len(tissue_pairs & subset_pairs)
        rows.append({
            "tissue": tissue,
            "n_total": len(tissue_pairs),
            "n_recovered": n_recovered,
            "n_missed": len(tissue_pairs) - n_recovered,
            "frac_recovered": n_recovered / len(tissue_pairs) if tissue_pairs else float("nan"),
        })
    summary = pd.DataFrame(rows).sort_values("n_total", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7, 0.6 * len(summary) + 1.5))
    y = range(len(summary))
    bar_height = 0.55

    ax.barh(y, summary["n_recovered"], height=bar_height, color=COLOR_RECOVERED,
            label="Also an outlier in Subset", zorder=3)
    ax.barh(y, summary["n_missed"], height=bar_height, left=summary["n_recovered"],
            color=COLOR_MISSED, label="Not an outlier in Subset", zorder=3)

    # direct labels: recovered fraction, placed just past the end of each bar
    for yi, row in zip(y, summary.itertuples()):
        ax.text(row.n_total + row.n_total * 0.015, yi, f"{row.frac_recovered:.0%} ({row.n_recovered}/{row.n_total})",
                va="center", ha="left", fontsize=9, color=COLOR_TEXT)

    ax.set_yticks(list(y))
    ax.set_yticklabels(summary["tissue"])
    ax.set_xlabel("Sample–junction outlier pairs")
    ax.set_title("Tissue outliers also flagged in the joint run")

    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.xaxis.grid(True, color=COLOR_GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(left=False)

    ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()

    return fig, ax, summary
