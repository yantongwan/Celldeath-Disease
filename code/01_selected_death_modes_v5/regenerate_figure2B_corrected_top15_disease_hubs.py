#!/usr/bin/env python3
"""Regenerate corrected Figure 2B with the top 15 original disease hubs.

This standalone panel uses the existing selected death-mode 2000-2025 hub source table.
It is intended for replacing the Figure 2B panel in Illustrator when the panel
layout displays fewer than the intended 15 disease terms.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_figure2b_corrected")

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("CELLDEATH_FIGURE_ROOT", Path(__file__).resolve().parents[2]))
TABLE_DIR = ROOT / "table" / "02_Figure_2_breadth_hubs_journal_country_radar"
OUT_DIR = ROOT / "corrected_panels" / "Figure2B_top15_original_disease_hubs"
SOURCE = TABLE_DIR / "hub_score_original_selected_death_modes_2000_2025.csv"

STEM = "Figure2B_corrected_top15_original_disease_hubs"

ONCOLOGY_COLOR = "#B9553C"
NON_ONCOLOGY_COLOR = "#6C9EAE"
TEXT_COLOR = "#263238"
GRID_COLOR = "#E7E7E7"


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 160,
        "savefig.dpi": 600,
        "font.size": 7,
        "axes.titlesize": 8.4,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)


def compact_label(value: str, width: int = 34) -> str:
    value = str(value)
    if len(value) <= width:
        return value
    return value[: width - 3] + "..."


def load_top15() -> pd.DataFrame:
    if not SOURCE.exists():
        raise FileNotFoundError(f"Missing source table: {SOURCE}")
    df = pd.read_csv(SOURCE)
    required = {
        "hub_rank",
        "disease_term",
        "disease_system",
        "death_mode_breadth",
        "pair_count",
        "oncology_flag",
        "hub_score",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Hub source table is missing columns: {sorted(missing)}")
    top15 = df.sort_values("hub_score", ascending=False).head(15).copy()
    top15["display_rank"] = np.arange(1, len(top15) + 1)
    if len(top15) != 15:
        raise ValueError(f"Expected 15 rows for Figure 2B, found {len(top15)}")
    return top15


def save_all(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(
        OUT_DIR / f"{stem}.tiff",
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def draw_panel(top15: pd.DataFrame) -> None:
    plot_df = top15.sort_values("hub_score", ascending=True).reset_index(drop=True)
    y = np.arange(len(plot_df))
    colors = [ONCOLOGY_COLOR if bool(v) else NON_ONCOLOGY_COLOR for v in plot_df["oncology_flag"]]

    fig, ax = plt.subplots(figsize=(4.55, 4.95))
    ax.barh(y, plot_df["hub_score"], color=colors, height=0.74, edgecolor="white", linewidth=0.35)

    labels = [
        f"{int(rank)}. {compact_label(term)}"
        for rank, term in zip(plot_df["display_rank"], plot_df["disease_term"])
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Literature hub score")
    ax.set_title("B  Original disease-literature hubs", loc="left", weight="bold")
    ax.grid(axis="x", color=GRID_COLOR, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(11.6, float(plot_df["hub_score"].max()) * 1.16))

    xmax = ax.get_xlim()[1]
    for yy, score, breadth, pairs in zip(
        y,
        plot_df["hub_score"],
        plot_df["death_mode_breadth"],
        plot_df["pair_count"],
    ):
        ax.text(
            score + xmax * 0.012,
            yy,
            f"{score:.2f} | {int(breadth)} modes, {int(pairs):,} pairs",
            va="center",
            ha="left",
            fontsize=5.45,
            color=TEXT_COLOR,
        )

    handles = [
        mpatches.Patch(color=ONCOLOGY_COLOR, label="Neoplasm/oncology term"),
        mpatches.Patch(color=NON_ONCOLOGY_COLOR, label="Other disease term"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.10),
        ncol=2,
        fontsize=5.8,
        handlelength=1.0,
        columnspacing=1.2,
    )
    ax.text(
        0,
        -0.24,
        "Top 15 by hub score. Hub score = z[log1p(pair count)] + z[death-mode breadth] + z[weighted MCS] + z[weighted DSS].",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.4,
        color="#555555",
    )
    fig.subplots_adjust(left=0.42, right=0.965, top=0.93, bottom=0.25)
    save_all(fig, STEM)


def write_source_and_qc(top15: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_source = OUT_DIR / f"{STEM}_source_data.csv"
    top15.to_csv(out_source, index=False)
    table_copy = TABLE_DIR / f"{STEM}_source_data.csv"
    top15.to_csv(table_copy, index=False)

    qc = pd.DataFrame(
        [
            {"metric": "source_table", "value": str(SOURCE)},
            {"metric": "intended_panel", "value": "Figure 2B"},
            {"metric": "display_rule", "value": "top 15 rows sorted by hub_score descending"},
            {"metric": "row_count", "value": len(top15)},
            {"metric": "first_disease", "value": top15.iloc[0]["disease_term"]},
            {"metric": "last_disease", "value": top15.iloc[-1]["disease_term"]},
        ]
    )
    qc.to_csv(OUT_DIR / f"{STEM}_QC.csv", index=False)
    qc.to_csv(TABLE_DIR / f"{STEM}_QC.csv", index=False)


def main() -> int:
    top15 = load_top15()
    write_source_and_qc(top15)
    draw_panel(top15)
    print(f"Wrote corrected Figure 2B panel to: {OUT_DIR}")
    print(f"Rows in corrected panel: {len(top15)}")
    print("Diseases:")
    for row in top15.itertuples(index=False):
        print(f"{int(row.display_rank):02d}. {row.disease_term}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
