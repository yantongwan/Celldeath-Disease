#!/usr/bin/env python3
"""Draw top journal and country death-mode stacked bars by PMID number.

This is the count-scale companion to the proportion figure. Segment colours
match Figure 2, while segment widths represent PMID/atlas counts.
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_selected_death_modes")

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
SUPP_SRC = V5 / "04_supplementary_figures_regenerated" / "source_data"
MAIN = V5 / "03_main_figures_regenerated"
STEM = "Figure_top_journal_country_death_mode_pmid_numbers_selected_death_modes_2000_2025"

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

PLOT_N = 15

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


def compact_label(value: str, width: int) -> str:
    value = str(value)
    if len(value) <= width:
        return value
    return textwrap.shorten(value, width=width, placeholder="...")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.075, 1.04, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))


def read_count_table(path: Path, entity_col: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing source table: {path}")
    df = pd.read_csv(path)
    missing = [c for c in [entity_col, *DEATH_ORDER] if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} missing columns: {missing}")
    for mode in DEATH_ORDER:
        df[mode] = pd.to_numeric(df[mode], errors="coerce").fillna(0).astype(float)
    df["atlas_count"] = df[DEATH_ORDER].sum(axis=1)
    df = df[df["atlas_count"] > 0].copy()
    df = df.sort_values("atlas_count", ascending=False).head(20)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


def make_long(df: pd.DataFrame, entity_col: str, panel: str) -> pd.DataFrame:
    out = df.melt(
        id_vars=[entity_col, "atlas_count", "rank"],
        value_vars=DEATH_ORDER,
        var_name="death_mode",
        value_name="death_mode_pmid_count",
    )
    out["death_mode_proportion"] = np.where(
        out["atlas_count"] > 0,
        out["death_mode_pmid_count"] / out["atlas_count"],
        0.0,
    )
    out["panel"] = panel
    out = out.rename(columns={entity_col: "entity"})
    return out


def plot_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    entity_col: str,
    title: str,
    label: str,
    wrap_width: int,
    plot_n: int = PLOT_N,
) -> None:
    plot_df = df.sort_values("atlas_count", ascending=False).head(plot_n)
    plot_df = plot_df.sort_values("atlas_count", ascending=True).reset_index(drop=True)
    y = np.arange(len(plot_df))
    left = np.zeros(len(plot_df), dtype=float)

    for mode in DEATH_ORDER:
        width = plot_df[mode].to_numpy(dtype=float)
        ax.barh(
            y,
            width,
            left=left,
            color=MODE_COLORS[mode],
            edgecolor="white",
            linewidth=0.25,
            height=0.78,
            label=mode,
        )
        left += width

    labels = [compact_label(v, wrap_width) for v in plot_df[entity_col]]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    xmax = float(plot_df["atlas_count"].max())
    ax.set_xlim(0, xmax * 1.15)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))
    for yy, total in zip(y, plot_df["atlas_count"]):
        ax.text(total + xmax * 0.018, yy, f"n={int(total):,}", va="center", ha="left", fontsize=5.8, color="#263238")
    ax.set_xlabel("PMID Number")
    ax.set_title(title, loc="left", weight="bold")
    panel_label(ax, label)
    clean_axes(ax)


def savefig(fig: plt.Figure, stem: str) -> None:
    for fmt in ["pdf", "png", "tiff"]:
        out = MAIN / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def main() -> int:
    journal = read_count_table(SUPP_SRC / "additional_journal_death_mode_counts_selected_death_modes_2000_2025.csv", "Journal")
    country = read_count_table(SUPP_SRC / "additional_country_death_mode_counts_selected_death_modes_2000_2025.csv", "Country")

    long = pd.concat(
        [
            make_long(journal, "Journal", "A Top journals by atlas count"),
            make_long(country, "Country", "B Country-level affiliation signal"),
        ],
        ignore_index=True,
    )

    source_dir = MAIN / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(source_dir / f"{STEM}_long_source_data.csv", index=False)
    journal.to_csv(source_dir / f"{STEM}_journal_wide_source_data.csv", index=False)
    country.to_csv(source_dir / f"{STEM}_country_wide_source_data.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 8.4), constrained_layout=False)
    plot_panel(axes[0], journal, "Journal", "Top journals by atlas count", "A", wrap_width=42)
    plot_panel(axes[1], country, "Country", "Country-level affiliation signal", "B", wrap_width=32)

    handles = [plt.Rectangle((0, 0), 1, 1, color=MODE_COLORS[m]) for m in DEATH_ORDER]
    fig.legend(
        handles,
        DEATH_ORDER,
        ncol=1,
        loc="center left",
        bbox_to_anchor=(0.835, 0.50),
        frameon=False,
        title="Cell death mode",
    )
    fig.suptitle(
        "Journal and country context of the eight-concept atlas",
        x=0.01,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    fig.text(
        0.01,
        0.020,
        "Bars are ranked by total atlas count; segment widths show PMID numbers by death mode. Country signal is affiliation-derived and should not be interpreted as national performance.",
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#555555",
    )
    fig.subplots_adjust(left=0.34, right=0.81, top=0.91, bottom=0.08, hspace=0.38)
    savefig(fig, STEM)

    qc = pd.DataFrame(
        [
            {"metric": "top_journals_available_in_source", "value": len(journal)},
            {"metric": "top_countries_available_in_source", "value": len(country)},
            {"metric": "top_journals_plotted", "value": min(PLOT_N, len(journal))},
            {"metric": "top_countries_plotted", "value": min(PLOT_N, len(country))},
            {"metric": "journal_max_atlas_count", "value": int(journal["atlas_count"].max())},
            {"metric": "country_max_atlas_count", "value": int(country["atlas_count"].max())},
            {"metric": "x_axis", "value": "PMID Number"},
            {"metric": "death_mode_palette_source", "value": "Figure_2 MODE_COLORS"},
        ]
    )
    qc.to_csv(source_dir / f"{STEM}_QC.csv", index=False)
    print(f"Wrote {MAIN / 'pdf' / f'{STEM}.pdf'}")
    print(f"Wrote source data under {source_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
