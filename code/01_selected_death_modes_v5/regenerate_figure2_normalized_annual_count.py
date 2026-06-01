#!/usr/bin/env python3
"""Regenerate Figure 2 with normalized annual counts without overwriting originals.

Panel A normalizes each death mode's annual unique-PMID count to its own
mode-specific maximum. Panel B keeps the original v5 takeoff/peak-year summary.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_selected_death_modes")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
START_YEAR = 2000
END_YEAR = 2025

STEM = "Figure_2_final_selected_death_modes_2000_2025_normalized_annual_count"

DEATH_ORDER = [
    "Ferroptosis",
    "Pyroptosis",
    "NETosis",
    "Necroptosis",
    "Immunogenic cell death",
    "Cuproptosis",
    "PANoptosis",
    "Disulfidptosis",
]

MODE_COLORS = {
    "Ferroptosis": "#B9553C",
    "Pyroptosis": "#D28A2E",
    "NETosis": "#5C8F7B",
    "Necroptosis": "#6B6BAE",
    "Immunogenic cell death": "#A05B8F",
    "Cuproptosis": "#3E7CB1",
    "PANoptosis": "#1B9E77",
    "Disulfidptosis": "#9C6A3D",
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
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "legend.fontsize": 6.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)
sns.set_context("paper")


def clean_axes(ax, axis: str = "x") -> None:
    ax.grid(axis=axis, color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)


def panel_label(ax, label: str) -> None:
    ax.text(-0.06, 1.06, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def savefig(fig: plt.Figure, stem: str) -> None:
    base = V5 / "03_main_figures_regenerated"
    for fmt in ["pdf", "png", "tiff"]:
        out = base / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def read_source_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    source_dir = V5 / "03_main_figures_regenerated" / "source_data"
    yearly_path = source_dir / "death_yearly_selected_death_modes_2000_2025.csv"
    summary_path = source_dir / "temporal_summary_selected_death_modes_2000_2025.csv"
    if not yearly_path.exists():
        raise FileNotFoundError(f"Missing yearly source data: {yearly_path}")
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing temporal summary source data: {summary_path}")
    yearly = pd.read_csv(yearly_path)
    summary = pd.read_csv(summary_path)
    required_yearly = {"Death_Mode", "Year", "Publication_Count"}
    required_summary = {"Death_Mode", "Introduction_Year", "Takeoff_Year", "Peak_Year"}
    if missing := required_yearly - set(yearly.columns):
        raise ValueError(f"Yearly source data missing columns: {sorted(missing)}")
    if missing := required_summary - set(summary.columns):
        raise ValueError(f"Temporal summary source data missing columns: {sorted(missing)}")
    return yearly, summary


def add_normalized_counts(yearly: pd.DataFrame) -> pd.DataFrame:
    yearly = yearly.copy()
    yearly["Publication_Count"] = pd.to_numeric(yearly["Publication_Count"], errors="coerce").fillna(0)
    mode_max = yearly.groupby("Death_Mode")["Publication_Count"].transform("max")
    yearly["Normalized_Annual_Count"] = np.where(mode_max > 0, yearly["Publication_Count"] / mode_max, 0.0)
    return yearly


def figure2_normalized(yearly: pd.DataFrame, summary: pd.DataFrame) -> None:
    yearly = add_normalized_counts(yearly)
    source_out = V5 / "03_main_figures_regenerated" / "source_data" / f"{STEM}_source_data.csv"
    yearly.to_csv(source_out, index=False)

    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.25, 0.75])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    panel_label(ax0, "A")
    eras = [
        (2000, 2006, "early disease-facing\nRCD vocabulary"),
        (2007, 2014, "inflammatory death\nadoption"),
        (2015, 2020, "ferroptosis/redox\nexpansion"),
        (2021, 2025, "metal-stress and\nintegrative phase"),
    ]
    era_cols = ["#F4F1E8", "#F8EBDD", "#EEF4F8", "#EAF3EF"]
    for (x0, x1, label), col in zip(eras, era_cols):
        ax0.axvspan(x0, x1, color=col, zorder=0)
        ax0.text(
            (x0 + x1) / 2,
            1.08,
            label,
            ha="center",
            va="bottom",
            fontsize=6.0,
            color="#555555",
            transform=ax0.get_xaxis_transform(),
        )

    for mode in DEATH_ORDER:
        sub = yearly[yearly["Death_Mode"] == mode].sort_values("Year")
        ax0.plot(sub["Year"], sub["Normalized_Annual_Count"], lw=1.45, color=MODE_COLORS[mode], label=mode)

    ax0.set_xlim(START_YEAR, END_YEAR)
    ax0.set_ylim(0, 1.08)
    ax0.set_ylabel("Normalized annual count")
    ax0.set_title("Normalized annual trajectory in disease literature", loc="left", weight="bold")
    ax0.legend(ncol=3, bbox_to_anchor=(0, -0.05), loc="upper left", frameon=False)
    clean_axes(ax0)

    panel_label(ax1, "B")
    summary = summary.set_index("Death_Mode").reindex(DEATH_ORDER).reset_index()
    y = np.arange(len(summary))[::-1]
    ax1.hlines(y, summary["Introduction_Year"], summary["Peak_Year"], color="#D0D0D0", lw=2)
    ax1.scatter(
        summary["Introduction_Year"],
        y,
        s=28,
        color="#FFFFFF",
        edgecolor="#555555",
        zorder=3,
        label="first year in filtered atlas",
    )
    ax1.scatter(summary["Takeoff_Year"], y, s=38, color=[MODE_COLORS[m] for m in summary["Death_Mode"]], zorder=3, label="takeoff")
    ax1.scatter(summary["Peak_Year"], y, s=24, marker="s", color="#333333", zorder=3, label="peak")
    ax1.set_yticks(y)
    ax1.set_yticklabels(summary["Death_Mode"])
    ax1.set_xlim(START_YEAR, END_YEAR + 1)
    ax1.set_xlabel("Year")
    ax1.set_title("Takeoff years after 2000-2025 filtering", loc="left", weight="bold")
    ax1.legend(ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, borderaxespad=0)
    clean_axes(ax1)

    fig.suptitle(
        "Figure 2. Historical waves of RCD terminology adoption, 2000-2025",
        x=0.01,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    savefig(fig, STEM)


def main() -> int:
    yearly, summary = read_source_tables()
    figure2_normalized(yearly, summary)
    print(f"Wrote normalized Figure 2 variant: {V5 / '03_main_figures_regenerated' / 'pdf' / f'{STEM}.pdf'}")
    print(f"Original Figure 2 was not overwritten: {V5 / '03_main_figures_regenerated' / 'pdf' / 'Figure_2_final_selected_death_modes_2000_2025.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
