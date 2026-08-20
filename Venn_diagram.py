import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
try:
    from matplotlib_venn import venn2
    _HAS_VENN = True
except ImportError:
    _HAS_VENN = False

CHART_PALETTE = {
    "blue":    "#2a78d6",
    "green":   "#008300",
    "magenta": "#e87ba4",
    "yellow":  "#eda100",
    "aqua":    "#1baf7a",
    "orange":  "#eb6834",
    "violet":  "#4a3aa7",
    "red":     "#e34948",
}
CATEGORICAL_ORDER = ["blue", "green", "magenta", "yellow", "aqua", "orange", "violet", "red"]

CHART_INK = {
    "surface":   "#fcfcfb",
    "primary":   "#0b0b0b",   # titles, value labels
    "secondary": "#52514e",  # subtitles, axis titles
    "muted":     "#898781",  # tick labels
    "gridline":  "#e1e0d9",  # hairline gridlines
    "baseline":  "#c3c2b7",  # axis/baseline
}


def apply_chart_style():
    """Set matplotlib rcParams to this notebook's shared chart style. Call once per figure."""
    plt.rcParams.update({
        "text.color": CHART_INK["primary"],
        "axes.edgecolor": CHART_INK["baseline"],
        "axes.labelcolor": CHART_INK["secondary"],
        "xtick.color": CHART_INK["muted"],
        "ytick.color": CHART_INK["muted"],
        "axes.facecolor": CHART_INK["surface"],
        "figure.facecolor": CHART_INK["surface"],
        "font.size": 10,
    })


def style_bar_axes(ax):
    """Strip an axes down to the bare style used for this notebook's bar charts:
    no chartjunk box, no gridlines, no y-axis tick numbers (values are direct-labeled
    on the bars instead), no tick marks."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(CHART_INK["baseline"])
    ax.spines["bottom"].set_color(CHART_INK["baseline"])
    ax.set_yticks([])
    ax.tick_params(axis="both", length=0)

def _load(src):
    """Accept either a path to a CSV/TSV or an already-loaded DataFrame."""
    if isinstance(src, pd.DataFrame):
        df = src.copy()
    else:
        ext = os.path.splitext(str(src))[1].lower()
        sep = '\t' if ext in ('.tsv', '.tab', '.txt') else ','
        df = pd.read_csv(src, sep=sep)
    df.rename(columns={'hgnc_symbol': 'hgncSymbol',
                       'JUNCTION_PADJ': 'padjust'}, inplace=True)
    return df


def compare_outputs(source_1, source_2, level='junction', cohort='',
                    name_1='R', name_2='Python', title_suffix='', out_path=None):
    """Compare the gene/junction sets of two outputs.

    source_1 / source_2 : path to a CSV or a DataFrame
    name_1 / name_2     : labels used in the plot and the printouts
    out_path             : if given, the figure is saved to this path
    """
    df_1 = _load(source_1)
    df_2 = _load(source_2)

    if level == 'gene':
        KEY_COLS = ['sampleID', 'hgncSymbol']
    elif level == 'junction':
        KEY_COLS = ['seqnames', 'start', 'end', 'sampleID', 'strand']
    else:
        raise ValueError(f"level must be 'gene' or 'junction', got {level!r}")

    keys_1 = set(df_1[KEY_COLS].astype(str).agg('|'.join, axis=1))
    keys_2 = set(df_2[KEY_COLS].astype(str).agg('|'.join, axis=1))

    only_1 = keys_1 - keys_2
    only_2 = keys_2 - keys_1
    shared = keys_1 & keys_2
    n_only_1, n_only_2, n_shared = len(only_1), len(only_2), len(shared)

    # Fixed categorical slots (position, not magnitude, decides the color):
    # left lobe = only source 1, lens = shared, right lobe = only source 2.
    c_1 = CHART_PALETTE['blue']
    c_shared = CHART_PALETTE['green']
    c_2 = CHART_PALETTE['magenta']

    apply_chart_style()
    fig, ax = plt.subplots(figsize=(6, 5))

    if _HAS_VENN:
        # venn2 subset order is (only A, only B, A∩B)
        v = venn2(subsets=(n_only_1, n_only_2, n_shared),
                  set_labels=(name_1, name_2), ax=ax)

        for region_id, color, count in (('10', c_1, n_only_1),
                                        ('01', c_2, n_only_2),
                                        ('11', c_shared, n_shared)):
            patch = v.get_patch_by_id(region_id)
            if patch is not None:                 # empty regions are not drawn
                patch.set_facecolor(color)
                patch.set_edgecolor('none')
                patch.set_alpha(0.75)
            label = v.get_label_by_id(region_id)
            if label is not None:
                label.set_text(f'{count:,}')
                label.set_fontsize(11)
                label.set_color(CHART_INK['primary'])

        for set_id in ('A', 'B'):
            label = v.get_label_by_id(set_id)
            if label is not None:
                label.set_fontsize(11)
                label.set_color(CHART_INK['secondary'])
    else:
        # Fallback: equally sized circles, overlap shown by alpha blending.
        r, dx = 1.0, 0.62
        ax.add_patch(Circle((-dx, 0), r, facecolor=c_1,
                            alpha=0.55, edgecolor='none', zorder=3))
        ax.add_patch(Circle((dx, 0), r, facecolor=c_2,
                            alpha=0.55, edgecolor='none', zorder=3))

        for x, count in ((-dx - r * 0.45, n_only_1),
                         (0, n_shared),
                         (dx + r * 0.45, n_only_2)):
            ax.text(x, 0, f'{count:,}', ha='center', va='center',
                    fontsize=11, color=CHART_INK['primary'], zorder=4)

        for x, label in ((-dx, name_1), (dx, name_2)):
            ax.text(x, r * 1.06, label, ha='center', va='bottom',
                    fontsize=11, color=CHART_INK['secondary'])

        ax.set_xlim(-dx - r * 1.2, dx + r * 1.2)
        ax.set_ylim(-r * 1.15, r * 1.4)
        ax.set_aspect('equal')

    ax.set_axis_off()

    total = max(len(keys_1), len(keys_2))
    pct_shared = 100 * n_shared / total if total else 0
    title = f'{level.capitalize()} set overlap'
    if cohort:
        title += f' in {cohort}'
    if title_suffix:
        title += f' {title_suffix}'
    ax.set_title(title, color=CHART_INK['primary'],
                 fontsize=13, fontweight='bold', loc='left', pad=26)
    ax.text(0, 1.03,
            f'{name_1}={len(keys_1):,}   {name_2}={len(keys_2):,}   '
            f'Shared={n_shared:,} ({pct_shared:.0f}%)',
            transform=ax.transAxes, ha='left', va='bottom', fontsize=9.5,
            color=CHART_INK['secondary'])

    plt.tight_layout()

    if out_path:
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches='tight')

    plt.show()

    return
