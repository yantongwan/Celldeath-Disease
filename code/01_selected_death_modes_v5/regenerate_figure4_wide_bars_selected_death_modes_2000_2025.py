#!/usr/bin/env python3
"""Regenerate Figure 4 with full-width hub-score bars.

The original v5 Figure 4 placed the original and no-neoplasm hub rankings
side-by-side. Long disease labels left each bar panel narrow, making rank
differences visually hard to inspect. This variant keeps the same data and
interpretation, but stacks panels vertically and annotates hub scores.
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_selected_death_modes")

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
MAIN = V5 / "03_main_figures_regenerated"
SRC = MAIN / "source_data"
STEM = "Figure_4_final_selected_death_modes_2000_2025_wide_bars"
ZOOM_STEM = "Figure_4_final_selected_death_modes_2000_2025_wide_bars_xstart5"
NO_LUNG_INJURY_ZOOM_STEM = "Figure_4_final_selected_death_modes_2000_2025_wide_bars_xstart5_no_lung_injury"
NO_LUNG_INJURY_TWO_COL_XSTART7_STEM = "Figure_4_final_selected_death_modes_2000_2025_two_col_xstart7_no_lung_injury"

SYSTEM_COLORS = {
    "Neoplasms": "#B9553C",
    "Metabolic and nutritional": "#A06A3B",
    "Respiratory": "#3E7CB1",
    "Musculoskeletal": "#8C6D62",
    "Pathologic process or sign": "#6FA4B3",
    "Cardiovascular": "#5C8F7B",
    "Endocrine": "#7A74B8",
    "Urogenital": "#6FA4B3",
    "Digestive": "#A05B8F",
    "Nervous system": "#6B6BAE",
    "Immune and inflammatory": "#D28A2E",
    "Other or mixed": "#8A8A8A",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 160,
        "savefig.dpi": 600,
        "font.size": 7,
        "axes.titlesize": 8.4,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)
sns.set_context("paper")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.055, 1.04, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def compact_label(value: str, width: int = 44) -> str:
    value = str(value)
    if len(value) <= width:
        return value
    return textwrap.shorten(value, width=width, placeholder="...")


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))


def read_source(name: str) -> pd.DataFrame:
    path = SRC / name
    if not path.exists():
        raise FileNotFoundError(f"Missing source data: {path}")
    return pd.read_csv(path)


def plot_hub_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    score_col: str,
    title: str,
    label: str,
    n: int,
    color_mode: str,
    x_min: float = 0.0,
    exclude_terms: set[str] | None = None,
) -> pd.DataFrame:
    if exclude_terms:
        frame = frame[~frame["disease_term"].isin(exclude_terms)].copy()
    plot_df = frame.sort_values(score_col, ascending=False).head(n).copy()
    plot_df = plot_df.sort_values(score_col, ascending=True).reset_index(drop=True)
    y = np.arange(len(plot_df))
    if color_mode == "oncology":
        colors = ["#B9553C" if bool(v) else "#6C9EAE" for v in plot_df["oncology_flag"]]
    else:
        colors = [SYSTEM_COLORS.get(v, "#8A8A8A") for v in plot_df["disease_system"]]

    ax.barh(y, plot_df[score_col], color=colors, height=0.82, edgecolor="white", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([compact_label(v) for v in plot_df["disease_term"]])
    xmax = max(12.0, float(plot_df[score_col].max()) * 1.13)
    ax.set_xlim(x_min, xmax)
    if x_min > 0:
        ax.axvline(x_min, color="#202020", lw=0.75)
    for yy, score, breadth, pairs in zip(y, plot_df[score_col], plot_df["death_mode_breadth"], plot_df["pair_count"]):
        ax.text(
            score + xmax * 0.010,
            yy,
            f"{score:.2f} | breadth={int(breadth)}, pairs={int(pairs):,}",
            va="center",
            ha="left",
            fontsize=5.7,
            color="#263238",
        )
    ax.set_xlabel("Literature hub score" + (f" (axis starts at {x_min:.1f})" if x_min > 0 else ""))
    ax.set_title(title, loc="left", weight="bold")
    panel_label(ax, label)
    clean_axes(ax)
    return plot_df


def add_metric_cards(ax: plt.Axes, oncology: pd.DataFrame, stats: pd.Series) -> None:
    ax.axis("off")
    cards = [
        (f"{int(oncology['oncology_flag'].sum()):,}", "oncology/neoplasm\nterms flagged", "#F7EAE6"),
        (f"{int((~oncology['oncology_flag']).sum()):,}", "non-oncology\nterms retained", "#EAF3EF"),
        (f"Spearman = {float(stats['spearman']):.3f}", "global ranking\nstructure", "#EEF4F8"),
        (f"Top25 Jaccard = {float(stats['top25_jaccard']):.3f}", "top hub identity\noverlap", "#F4F1E8"),
    ]
    for i, (big, small, color) in enumerate(cards):
        x = 0.02 + i * 0.245
        ax.add_patch(
            FancyBboxPatch(
                (x, 0.16),
                0.22,
                0.64,
                boxstyle="round,pad=0.012,rounding_size=0.02",
                fc=color,
                ec="#CCCCCC",
                lw=0.6,
            )
        )
        ax.text(x + 0.11, 0.56, big, ha="center", va="center", weight="bold", fontsize=8.6)
        ax.text(x + 0.11, 0.34, small, ha="center", va="center", fontsize=6.2)


def add_compact_metric_cards(ax: plt.Axes, oncology: pd.DataFrame, stats: pd.Series) -> None:
    ax.axis("off")
    cards = [
        (f"{int(oncology['oncology_flag'].sum()):,}", "oncology/neoplasm terms", "#F7EAE6"),
        (f"{int((~oncology['oncology_flag']).sum()):,}", "non-oncology terms", "#EAF3EF"),
        (f"Spearman {float(stats['spearman']):.3f}", "rank similarity", "#EEF4F8"),
        (f"Jaccard {float(stats['top25_jaccard']):.3f}", "top25 overlap", "#F4F1E8"),
    ]
    for i, (big, small, color) in enumerate(cards):
        x = 0.02 + i * 0.245
        ax.add_patch(
            FancyBboxPatch(
                (x, 0.20),
                0.22,
                0.56,
                boxstyle="round,pad=0.012,rounding_size=0.02",
                fc=color,
                ec="#CCCCCC",
                lw=0.55,
            )
        )
        ax.text(x + 0.11, 0.54, big, ha="center", va="center", weight="bold", fontsize=7.2)
        ax.text(x + 0.11, 0.35, small, ha="center", va="center", fontsize=5.5)


def savefig(fig: plt.Figure, stem: str) -> None:
    for fmt in ["pdf", "png", "tiff"]:
        out = MAIN / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def draw_figure(stem: str, x_min: float = 0.0, exclude_terms: set[str] | None = None) -> None:
    orig = read_source("hub_score_original_selected_death_modes_2000_2025.csv")
    no = read_source("hub_score_no_neoplasms_selected_death_modes_2000_2025.csv")
    oncology_path = V5 / "01_source_data_filtered" / "oncology_term_classification_selected_death_modes_2000_2025.csv"
    if not oncology_path.exists():
        raise FileNotFoundError(f"Missing oncology classification source: {oncology_path}")
    oncology = pd.read_csv(oncology_path)
    stats = read_source("oncology_rank_stats_selected_death_modes_2000_2025.csv").iloc[0]

    fig = plt.figure(figsize=(8.7, 9.3), constrained_layout=False)
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.35, 1.08, 0.24], hspace=0.34)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])

    orig_plot = plot_hub_panel(
        ax0,
        orig,
        "hub_score",
        "Original eight-concept hubs",
        "A",
        n=16,
        color_mode="oncology",
        x_min=x_min,
        exclude_terms=exclude_terms,
    )
    no_plot = plot_hub_panel(
        ax1,
        no,
        "hub_score_no_neoplasms",
        "Retained non-oncology structure after neoplasm removal",
        "B",
        n=12,
        color_mode="system",
        x_min=x_min,
        exclude_terms=exclude_terms,
    )
    add_metric_cards(ax2, oncology, stats)
    panel_label(ax2, "C")

    fig.suptitle(
        "Figure 4. Oncology-sensitive literature hubs in the eight-concept atlas",
        x=0.01,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    fig.text(
        0.01,
        0.012,
        "Hub score = z[log1p(pair count)] + z[death-mode breadth] + z[weighted cleaned MCS] + z[weighted DSS]. Scores describe literature structure, not biological convergence."
        + (f" A/B x-axes start at {x_min:.1f} to emphasize rank differences." if x_min > 0 else ""),
        ha="left",
        va="bottom",
        fontsize=5.9,
        color="#555555",
    )
    fig.subplots_adjust(left=0.36, right=0.98, top=0.93, bottom=0.07)
    savefig(fig, stem)

    out_dir = MAIN / "source_data"
    orig_plot.sort_values("hub_score", ascending=False).to_csv(out_dir / f"{stem}_panel_A_source_data.csv", index=False)
    no_plot.sort_values("hub_score_no_neoplasms", ascending=False).to_csv(out_dir / f"{stem}_panel_B_source_data.csv", index=False)
    pd.DataFrame(
        [
            {"component": "pair_count", "transformation": "log1p, then z-score", "interpretation": "pair-volume support"},
            {"component": "death_mode_breadth", "transformation": "z-score", "interpretation": "number of death modes linked to the disease term"},
            {"component": "mean_cleaned_mcs_proxy", "transformation": "Pair_Count-weighted mean, then z-score", "interpretation": "cleaned mechanism-vocabulary convergence proxy"},
            {"component": "mean_dss", "transformation": "Pair_Count-weighted mean, then z-score", "interpretation": "disease-literature specificity score"},
        ]
    ).to_csv(out_dir / f"{stem}_hub_score_formula.csv", index=False)
    print(f"Wrote {MAIN / 'pdf' / f'{stem}.pdf'}")
    print(f"Wrote source data under {out_dir}")


def draw_two_column_figure(stem: str, x_min: float = 7.0, exclude_terms: set[str] | None = None) -> None:
    orig = read_source("hub_score_original_selected_death_modes_2000_2025.csv")
    no = read_source("hub_score_no_neoplasms_selected_death_modes_2000_2025.csv")
    oncology_path = V5 / "01_source_data_filtered" / "oncology_term_classification_selected_death_modes_2000_2025.csv"
    if not oncology_path.exists():
        raise FileNotFoundError(f"Missing oncology classification source: {oncology_path}")
    oncology = pd.read_csv(oncology_path)
    stats = read_source("oncology_rank_stats_selected_death_modes_2000_2025.csv").iloc[0]

    fig = plt.figure(figsize=(8.2, 4.9), constrained_layout=False)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.20], width_ratios=[1, 1], hspace=0.32, wspace=0.62)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])

    orig_plot = plot_hub_panel(
        ax0,
        orig,
        "hub_score",
        "Original hubs",
        "A",
        n=14,
        color_mode="oncology",
        x_min=x_min,
        exclude_terms=exclude_terms,
    )
    no_plot = plot_hub_panel(
        ax1,
        no,
        "hub_score_no_neoplasms",
        "After neoplasm removal",
        "B",
        n=10,
        color_mode="system",
        x_min=x_min,
        exclude_terms=exclude_terms,
    )
    add_compact_metric_cards(ax2, oncology, stats)
    panel_label(ax2, "C")

    for ax in [ax0, ax1]:
        ax.tick_params(axis="y", labelsize=5.5)
        ax.tick_params(axis="x", labelsize=5.6)
        ax.xaxis.label.set_size(5.8)
        ax.title.set_size(7.0)
        for text in ax.texts:
            text.set_fontsize(4.7)

    fig.suptitle(
        "Figure 4. Oncology-sensitive literature hubs in the eight-concept atlas",
        x=0.01,
        ha="left",
        fontsize=9,
        weight="bold",
    )
    fig.text(
        0.01,
        0.012,
        f"Display excludes Lung Injury; next-ranked terms fill the panels. A/B x-axes start at {x_min:.1f}. Hub scores describe literature structure, not biological convergence.",
        ha="left",
        va="bottom",
        fontsize=5.3,
        color="#555555",
    )
    fig.subplots_adjust(left=0.23, right=0.985, top=0.88, bottom=0.11)
    savefig(fig, stem)

    out_dir = MAIN / "source_data"
    orig_plot.sort_values("hub_score", ascending=False).to_csv(out_dir / f"{stem}_panel_A_source_data.csv", index=False)
    no_plot.sort_values("hub_score_no_neoplasms", ascending=False).to_csv(out_dir / f"{stem}_panel_B_source_data.csv", index=False)
    pd.DataFrame(
        [
            {"display_rule": "excluded disease_term", "value": "; ".join(sorted(exclude_terms or []))},
            {"display_rule": "x_axis_start", "value": x_min},
            {"display_rule": "panel_A_rows", "value": len(orig_plot)},
            {"display_rule": "panel_B_rows", "value": len(no_plot)},
            {"display_rule": "interpretation", "value": "literature-structure descriptor, not biological convergence"},
        ]
    ).to_csv(out_dir / f"{stem}_display_QC.csv", index=False)
    print(f"Wrote {MAIN / 'pdf' / f'{stem}.pdf'}")
    print(f"Wrote source data under {out_dir}")


def main() -> int:
    draw_figure(STEM, x_min=0.0)
    draw_figure(ZOOM_STEM, x_min=5.0)
    draw_figure(NO_LUNG_INJURY_ZOOM_STEM, x_min=5.0, exclude_terms={"Lung Injury"})
    draw_two_column_figure(NO_LUNG_INJURY_TWO_COL_XSTART7_STEM, x_min=7.0, exclude_terms={"Lung Injury"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
