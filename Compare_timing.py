import csv
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE = Path("/s/project/py_fraser/gtex/output_gtex")

R_STEP_MAP = {
    "Initializing dataset (load + jaccard calculation)": "load+jaccard",
    "Filtering expression and variability": "filter",
    "Finding latent dimension": "find_q",
    "Fitting autoencoder": "fit_autoencoder",
    "Computing p-values": "compute_pvalues",
    "Annotating genes": "annotate_genes",
    "Adjusting p-values": "adjust_pvalues",
    "Extracting results": "extract_results",
}
PY_STEP_MAP = {
    "Initializing dataset": "load+jaccard",  # in Python this also includes filtering + GTF annotation
    "Finding latent dimension": "find_q",
    "Initializing model": "init_model",
    "Computing initial loss": "init_loss",
    "Computing final loss": "final_loss",
    "Computing p-values and residuals": "fit+compute_pvalues",  # older logs: fit + p-values were one step
    "Fitting residuals": "fit_residuals",  # newer logs: split into two steps
    "Computing p-values": "compute_pvalues",
    "Adjusting p-values": "adjust_pvalues",
    "Finalizing model": "finalize_model",
}

# Aligned buckets: same real-world phase, regardless of how each side times it.
# Keys inside a bucket are RUN KINDS, so every run of the same kind reuses them.
BUCKETS = {
    "load_filter_annotate": {"R": ["load+jaccard", "filter", "annotate_genes"], "Python": ["load+jaccard"]},
    "find_latent_dim": {"R": ["find_q"], "Python": ["find_q"]},
    "fit_autoencoder": {
        "R": ["fit_autoencoder"],
        # Older logs never split "fit+compute_pvalues" into its two parts, so that
        # combined time is counted entirely here, leaving compute_pvalues at 0 for those runs.
        "Python": ["init_model", "init_loss", "final_loss", "fit_residuals", "fit+compute_pvalues"],
    },
    "compute_pvalues": {"R": ["compute_pvalues"], "Python": ["compute_pvalues"]},
    "adjust_pvalues": {"R": ["adjust_pvalues"], "Python": ["adjust_pvalues"]},
    "extract_save_results": {"R": ["extract_results"], "Python": ["finalize_model", "save_output"]},
}

BUCKET_LABELS = {
    "load_filter_annotate": "Load + filter + annotate",
    "find_latent_dim": "Find latent dim",
    "fit_autoencoder": "Fit",
    "compute_pvalues": "Compute p-values",
    "adjust_pvalues": "Adjust p-values",
    "extract_save_results": "Extract / save results",
}

# Keys a parser may return that are metadata, not step durations.
META_KEYS = ("n_samples", "n_junctions_before", "n_junctions_after")

PY_LINE_RE = re.compile(r"^(?P<ts>\S+ \S+) - \S+ - INFO - (?P<step>.+?): (?P<secs>[\d.]+)s \(total")
PY_TS_RE = re.compile(r"^(?P<ts>\S+ \S+) - ")

# Junction/sample counts, logged once per run while building the FraserDataset.
PY_N_SAMPLES_RE = re.compile(r"Number of samples in the dataset: (?P<n>\d+)")
PY_EXPR_FILTER_RE = re.compile(r"Junctions passing expression filter: (?P<passed>\d+) / (?P<total>\d+)")
PY_VAR_FILTER_RE = re.compile(r"Junctions passing variability filter: (?P<passed>\d+) / (?P<total>\d+)")

DEFAULT_COLORS = ["#6D6D6D", "#00D223", "#004719", "#B65C00", "#1F5FA8", "#8E44AD", "#C0392B", "#00838F"]


@dataclass(frozen=True)
class Run:
    """One labelled result tree to compare.

    label:   shown in the legend / CSV / subtitle; must be unique.
    kind:    "R" or "Python" — selects the parser and the bucket step lists.
    pattern: path to one result file, with `{cohort}` marking the cohort part.
             Relative to `base_dir`, or absolute. Globs (`*`) are allowed
             elsewhere in the path, e.g. "R_V2/AE_CPU/{cohort}/*_timing.csv".
    color:   optional override; otherwise taken from DEFAULT_COLORS.
    """

    label: str
    kind: str
    pattern: str
    color: str = None


# ---------------------------------------------------------------- parsers


def _parse_r(csv_path):
    out = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            out[R_STEP_MAP.get(row["step"], row["step"])] = float(row["seconds"])
    return out


def _parse_py(log_path):
    # A job can be requeued/restarted with the same log file appended to rather
    # than truncated — anchor on the LAST "Starting protrider" so idle time
    # before a restart doesn't inflate the wall time.
    out = {}
    first_ts = last_ts = saving_start_ts = None
    with open(log_path) as f:
        for line in f:
            m = PY_TS_RE.match(line)
            if m:
                ts = datetime.strptime(m.group("ts"), "%Y-%m-%d %H:%M:%S,%f")
                if "Starting protrider" in line:
                    first_ts = ts
                    out.clear()
                    saving_start_ts = None
                elif first_ts is None:
                    first_ts = ts
                last_ts = ts
            if "Saving results in long format" in line:
                saving_start_ts = ts
            sm = PY_LINE_RE.match(line)
            if sm:
                out[PY_STEP_MAP.get(sm.group("step"), sm.group("step"))] = float(sm.group("secs"))

            # Sample count and junction counts before/after filtering.
            nm = PY_N_SAMPLES_RE.search(line)
            if nm:
                out["n_samples"] = int(nm.group("n"))
            em = PY_EXPR_FILTER_RE.search(line)
            if em:
                out["n_junctions_before"] = int(em.group("total"))
            vm = PY_VAR_FILTER_RE.search(line)
            if vm:
                out["n_junctions_after"] = int(vm.group("passed"))
    if saving_start_ts and last_ts:
        out["save_output"] = (last_ts - saving_start_ts).total_seconds()
    if first_ts and last_ts:
        out["TOTAL_WALLTIME"] = (last_ts - first_ts).total_seconds()
    return out


PARSERS = {"R": _parse_r, "Python": _parse_py}


# ---------------------------------------------------------------- discovery


def _discover_paths(base, run):
    """Map cohort -> path for one run, by globbing its pattern."""
    pattern = run.pattern
    if os.path.isabs(pattern):
        root, rel_pattern = Path("/"), pattern.lstrip("/")
    else:
        root, rel_pattern = Path(base), pattern

    rx = re.compile(
        re.escape(rel_pattern)
        .replace(r"\{cohort\}", r"(?P<cohort>[^/]+)")
        .replace(r"\*", r"[^/]*")
        + "$"
    )
    found = {}
    for p in sorted(root.glob(rel_pattern.format(cohort="*"))):
        m = rx.match(p.relative_to(root).as_posix())
        if m:
            found[m.group("cohort")] = p
    return found


def _total_seconds(parsed):
    if "TOTAL_WALLTIME" in parsed:
        return parsed["TOTAL_WALLTIME"]
    return sum(v for k, v in parsed.items() if k not in META_KEYS and k != "TOTAL_WALLTIME")


# ---------------------------------------------------------------- comparison


def compare_timing(runs, base_dir=BASE, cohorts=None, csv_out=None, require_all=True):
    """Compare aligned-bucket step timing across any number of runs.

    Returns LONG-format rows: one dict per (cohort, bucket, run) with keys
    {cohort, bucket, run, kind, seconds}. The bucket == "TOTAL" rows also carry
    n_junctions_before/after and n_samples when the run's parser found them.

    require_all=True keeps only cohorts present in every run (fair comparison);
    set it to False to keep any cohort found in at least one run.
    """
    base = Path(base_dir)
    paths = {r.label: _discover_paths(base, r) for r in runs}

    if cohorts:
        selected = sorted(cohorts)
    elif require_all:
        selected = sorted(set.intersection(*(set(p) for p in paths.values())) if paths else [])
    else:
        selected = sorted(set().union(*(set(p) for p in paths.values())) if paths else [])

    rows = []
    for cohort in selected:
        parsed = {}
        for run in runs:
            path = paths[run.label].get(cohort)
            if path is None:
                print(f"Skipping {cohort} / {run.label}: no file matching {run.pattern}")
                continue
            parsed[run.label] = PARSERS[run.kind](path)

        for bucket, srcs in BUCKETS.items():
            for run in runs:
                if run.label not in parsed:
                    continue
                d = parsed[run.label]
                rows.append({
                    "cohort": cohort,
                    "bucket": bucket,
                    "run": run.label,
                    "kind": run.kind,
                    "seconds": sum(d.get(s, 0.0) for s in srcs[run.kind]),
                })

        for run in runs:
            if run.label not in parsed:
                continue
            d = parsed[run.label]
            row = {
                "cohort": cohort,
                "bucket": "TOTAL",
                "run": run.label,
                "kind": run.kind,
                "seconds": _total_seconds(d),
            }
            row.update({k: d.get(k) for k in META_KEYS})
            rows.append(row)

    if csv_out:
        fields = ["cohort", "bucket", "run", "kind", "seconds", *META_KEYS]
        with open(csv_out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {csv_out}")

    return rows


def pivot_rows(rows, value="seconds", scale=1 / 60):
    """Long rows -> {(cohort, bucket): {run: value}}, handy for a wide table."""
    wide = {}
    for r in rows:
        wide.setdefault((r["cohort"], r["bucket"]), {})[r["run"]] = r[value] * scale
    return wide


# ---------------------------------------------------------------- plotting


def plot_timing(rows, runs=None, baseline=None, bucket_order=None, ncols=2,
                figsize_per_panel=(5.5, 3.2), title="FRASER step timing by cohort", output_path = ''):
    """Grouped horizontal bars of per-bucket timing, one panel per cohort, one bar per run.

    `runs` (the same list passed to compare_timing) fixes the legend order and
    colors; omit it and the order is taken from the rows. `baseline` is the run
    label speedups are computed against — defaults to the first run.
    """
    if runs:
        labels = [r.label for r in runs]
        colors = {r.label: (r.color or DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
                  for i, r in enumerate(runs)}
    else:
        labels = list(dict.fromkeys(r["run"] for r in rows))
        colors = {lab: DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i, lab in enumerate(labels)}
    baseline = baseline or labels[0]

    bucket_order = bucket_order or list(BUCKET_LABELS)
    cohorts = sorted({r["cohort"] for r in rows})
    nrows = -(-len(cohorts) // ncols)  # ceil

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    n = len(labels)
    height = 0.8 / n

    for ax, cohort in zip(axes, cohorts):
        vals = {(r["bucket"], r["run"]): r["seconds"] for r in rows if r["cohort"] == cohort}
        present = {b for b, _ in vals}
        buckets = [b for b in bucket_order if b in present]
        y = np.arange(len(buckets))

        for i, lab in enumerate(labels):
            offset = (i - (n - 1) / 2) * height  # first run ends up on top after inverting
            minutes = [vals.get((b, lab), 0.0) / 60 for b in buckets]
            ax.barh(y + offset, minutes, height=height, color=colors[lab], label=lab)

        ax.set_yticks(y)
        ax.set_yticklabels([BUCKET_LABELS.get(b, b) for b in buckets], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("minutes", fontsize=9, color="#52514e")
        ax.tick_params(axis="x", labelsize=8, colors="#52514e")
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#c3c2b7")
        ax.grid(axis="x", color="#e1e0d9", linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)

        totals = {r["run"]: r for r in rows if r["cohort"] == cohort and r["bucket"] == "TOTAL"}
        parts = []
        base_tot = totals.get(baseline, {}).get("seconds")
        for lab in labels:
            if lab not in totals:
                continue
            tot = totals[lab]["seconds"] / 60
            if base_tot and lab != baseline:
                parts.append(f"{lab} {tot:.0f}m ({base_tot / 60 / tot:.2f}x)")
            else:
                parts.append(f"{lab} {tot:.0f}m")
        # wrap two entries per line so the subtitle stays readable with many runs
        subtitle = "\n".join(" · ".join(parts[i:i + 2]) for i in range(0, len(parts), 2))

        meta = next((t for t in totals.values()
                     if t.get("n_junctions_before") and t.get("n_junctions_after")), None)
        if meta:
            subtitle += f"\njunctions {meta['n_junctions_before']:,} → {meta['n_junctions_after']:,}"
            if meta.get("n_samples") is not None:
                subtitle += f" · samples {meta['n_samples']:,}"
        ax.set_title(f"{cohort}\n{subtitle}", fontsize=10, loc="left")

    for ax in axes[len(cohorts):]:
        ax.axis("off")

    handles, plotted = axes[0].get_legend_handles_labels()
    fig.legend(handles, plotted, loc="upper right", frameon=False, fontsize=9)
    fig.suptitle("FRASER step timing by cohort", fontsize=13, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    if output_path:
        fig.savefig(output_path, dpi=300)
        print(f"Saved figure to {output_path}")

    return fig
    