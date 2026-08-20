import pandas as pd
import numpy as np
import plotnine as pn
from scipy.stats import beta
from datetime import date

USECOLS = {
    "python": ["sampleID", "seqnames", "start", "end", "strand", "JUNCTION_PVALUE", "JUNCTION_outlier"],
    "r": ["sampleID", "seqnames", "start", "end", "pValue"],  
}
READ_KWARGS = {
    "python": dict(compression="gzip"),
    "r": dict(compression="gzip", sep="\t"),
}
RENAME_COLS = {
    "python": {"JUNCTION_PVALUE": "pValue"},
    "r": {},
}


today = date.today().strftime("%Y-%m-%d")  # e.g. 2026-07-29

def add_variant_support(df, vep_path, sample_col="sampleID", chrom_col="seqnames",
                         start_col="start", end_col="end", window=500, impact="HIGH"):
    """
    Adds a boolean 'has_variant' column: True where the junction has a
    supporting variant — same subject, same chromosome, and variant position
    within [start - window, end + window].
    """
    vep = pd.read_csv(vep_path)
    vep = vep[vep["IMPACT"] == impact]

    variant_pos = vep["variantID"].str.split("_", expand=True)
    vep = vep.assign(chrom=variant_pos[0], pos=variant_pos[1].astype(int))

    df = df.reset_index(drop=True)
    subject_id = df[sample_col].str.split("-").str[:2].str.join("-")

    candidates = pd.DataFrame({
        "_row_id": df.index,
        "_subject_id": subject_id,
        "_chrom": df[chrom_col],
        "_start": df[start_col],
        "_end": df[end_col],
    }).merge(vep[["subjectID", "chrom", "pos"]],
             left_on=["_subject_id", "_chrom"], right_on=["subjectID", "chrom"])

    supporting = candidates[
        (candidates["pos"] >= candidates["_start"] - window) &
        (candidates["pos"] <= candidates["_end"] + window)
    ]

    df = df.copy()
    df["has_variant"] = df.index.isin(set(supporting["_row_id"]))
    return df


def calculate_splicing_proportions(df, rank_col="pValue", variant_col="has_variant"):
    """Sorts by rank_col (ascending) and computes, at each rank, the running
    proportion of outliers with a supporting variant plus a 95% CI."""
    df = df.sort_values(rank_col, ascending=True).reset_index(drop=True)

    is_hit = df[variant_col].to_numpy()
    rank = np.arange(1, len(df) + 1)
    cum_hits = np.cumsum(is_hit)

    alpha = 0.05
    df["rank"] = rank
    df["proportion"] = cum_hits / rank
    df["ci_min"] = beta.ppf(alpha / 2, cum_hits, rank - cum_hits + 1)
    df["ci_max"] = beta.ppf(1 - alpha / 2, cum_hits + 1, rank - cum_hits)
    return df



def plot_splicing_enrichment(df, rank_col="rank", prop_col="proportion",
                              ci_min_col="ci_min", ci_max_col="ci_max",
                              color_col="source", min_rank=500, max_rank=1000000,
                              show_outlier_count = True, outlier_source = "Python",
                              outlier_col="JUNCTION_outlier",
                              ylabel="Proportion of splicing outliers with\nrare splice-disrupting variants",
                              out_path=None):
    plot_data = df[(df[rank_col] > min_rank) & (df[rank_col] < max_rank)]

    p = (
        pn.ggplot(plot_data)
        + pn.geom_line(pn.aes(x=rank_col, y=prop_col, color=color_col))
        + pn.geom_ribbon(pn.aes(x=rank_col, ymin=ci_min_col, ymax=ci_max_col, fill=color_col), alpha=0.2, outline_type='none')
        + pn.scale_x_log10()
        + pn.annotation_logticks(sides="b")
        + pn.labs(y=ylabel, x="Outlier rank cutoff")
        + pn.theme_bw(base_size=12)
    )

    if show_outlier_count:
            outlier_df = df if outlier_source is None else df[df[color_col] == outlier_source]
            outlier_count = int(outlier_df[outlier_col].sum())
            p = p + pn.geom_vline(xintercept=outlier_count, linetype="dashed", color="black")

    
    if out_path:
        p.save(out_path, width=6, height=4, dpi=300, bbox_inches='tight')
    return p

def load_junctions(path, kind):
    """Reads a junctions table from either a Python (kind='python') or R
    (kind='r') FRASER output file and normalizes it to a common schema."""
    df = pd.read_csv(path, usecols=USECOLS[kind], **READ_KWARGS[kind])
    return df.rename(columns=RENAME_COLS[kind])


def prepare_junctions(df, vep_path, label, rank_col="pValue"):
    df = add_variant_support(df, vep_path=vep_path)
    df = calculate_splicing_proportions(df, rank_col=rank_col)
    df["source"] = label
    print(f"Prepared {label} junctions: {len(df)} ({df['has_variant'].sum()} with variants)")
    return df


def compare_splicing_runs(runs, vep_path, out_path, ylabel, outlier_source,
                           rank_col="pValue"):
    """Compares an arbitrary set of FRASER runs (any mix of Python/R) on the
    same cohort and plots splicing-outlier enrichment for supporting variants.

    runs: list of dicts, each {"path": ..., "kind": "python"|"r", "label": ...}
    outlier_source: label (matching one run's "label") used for the outlier-count vline
    """
    frames = []
    for run in runs:
        df = load_junctions(run["path"], run["kind"])
        df = prepare_junctions(df, vep_path, run["label"], rank_col=rank_col)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    return plot_splicing_enrichment(combined, ylabel=ylabel,
                                     outlier_source=outlier_source, out_path=out_path)

