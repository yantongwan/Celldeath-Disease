#!/usr/bin/env python3
"""Redraw balanced evidence-stage signal map as disease radar profiles.

The original balanced dot plot is preserved. This script creates a less crowded
radar small-multiple view using the same balanced top-pair source table.
"""

from __future__ import annotations

import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_evidence_radar")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
SUPP = V5 / "04_supplementary_figures_regenerated"
SRC = SUPP / "source_data"

SOURCE = SRC / "additional_evidence_stage_signal_balanced_top_pairs_selected_death_modes_2000_2025.csv"
STEM = "Additional_Figure_evidence_stage_signal_radar_balanced_selected_death_modes_2000_2025"

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
        "xtick.labelsize": 5.2,
        "ytick.labelsize": 5.0,
        "legend.fontsize": 6.0,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)
sns.set_context("paper")


def compact_label(value: str, width: int = 20) -> str:
    value = str(value)
    if len(value) <= width:
        return value
    return textwrap.shorten(value, width=width, placeholder="...")


def wrapped_label(value: str, width: int = 16) -> str:
    value = compact_label(value, width=32)
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False))


def savefig(fig: plt.Figure, stem: str) -> None:
    for fmt in ["pdf", "png", "tiff"]:
        out = SUPP / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def load_balanced() -> pd.DataFrame:
    if not SOURCE.exists():
        raise FileNotFoundError(f"Missing source table: {SOURCE}")
    df = pd.read_csv(SOURCE)
    required = {
        "death_mode",
        "disease_term",
        "pair_count",
        "audited_stage_score_max",
        "clinical_vocabulary_fraction",
        "interventional_supported_article_n",
        "evidence_stage_signal",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"{SOURCE.name} missing columns: {missing}")
    df = df[df["death_mode"].isin(DEATH_ORDER)].copy()
    df["evidence_stage_signal"] = pd.to_numeric(df["evidence_stage_signal"], errors="coerce")
    df["pair_count"] = pd.to_numeric(df["pair_count"], errors="coerce")
    df = df.dropna(subset=["evidence_stage_signal", "pair_count"])
    df["death_mode"] = pd.Categorical(df["death_mode"], categories=DEATH_ORDER, ordered=True)
    df = df.sort_values(["death_mode", "evidence_stage_signal"], ascending=[True, False])
    df["within_mode_rank"] = df.groupby("death_mode", observed=False).cumcount() + 1
    return df


def plot_radar_panel(ax: plt.Axes, sub: pd.DataFrame, mode: str, rmax: float) -> None:
    sub = sub.sort_values("evidence_stage_signal", ascending=False).head(5).copy()
    labels = [wrapped_label(x) for x in sub["disease_term"]]
    values = sub["evidence_stage_signal"].to_numpy(dtype=float)
    n = len(values)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    values_closed = np.r_[values, values[0]]
    angles_closed = np.r_[angles, angles[0]]
    color = MODE_COLORS.get(mode, "#888888")

    ax.plot(angles_closed, values_closed, color=color, lw=1.7)
    ax.fill(angles_closed, values_closed, color=color, alpha=0.18)
    ax.scatter(angles, values, s=np.clip(np.sqrt(sub["pair_count"].to_numpy(dtype=float)) * 8, 18, 90), color=color, edgecolor="white", linewidth=0.4, zorder=3)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=5.2)
    ax.tick_params(axis="x", pad=2)
    ax.set_ylim(0, rmax)
    ax.set_yticks([0.5, 1.0, 1.5])
    ax.set_yticklabels(["0.5", "1.0", "1.5"], fontsize=4.7, color="#666666")
    ax.yaxis.grid(True, color="#D9D9D9", lw=0.45)
    ax.xaxis.grid(True, color="#EEEEEE", lw=0.35)
    ax.spines["polar"].set_color("#D0D0D0")
    ax.spines["polar"].set_linewidth(0.55)

    # Use a color swatch plus number only; disease labels do not repeat death-mode text.
    ax.text(
        0.02,
        1.05,
        f"Top {len(sub)}",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=6.0,
        weight="bold",
        color="#222222",
    )
    ax.add_patch(
        plt.Rectangle(
            (0.78, 1.035),
            0.14,
            0.035,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            clip_on=False,
        )
    )


def main() -> int:
    df = load_balanced()
    out_source = SRC / f"{STEM}_source_data.csv"
    df.to_csv(out_source, index=False)

    rmax = max(1.6, float(np.ceil(df["evidence_stage_signal"].max() * 10) / 10))

    fig, axes = plt.subplots(
        2,
        4,
        figsize=(12.8, 6.8),
        subplot_kw={"projection": "polar"},
        constrained_layout=False,
    )
    for ax, mode in zip(axes.flat, DEATH_ORDER):
        plot_radar_panel(ax, df[df["death_mode"].eq(mode)], mode, rmax)

    handles = [plt.Line2D([0], [0], color=MODE_COLORS[m], lw=5) for m in DEATH_ORDER]
    fig.legend(
        handles,
        DEATH_ORDER,
        ncol=1,
        loc="center left",
        bbox_to_anchor=(0.835, 0.50),
        frameon=False,
        title="Cell death mode color",
    )
    fig.suptitle(
        "Balanced evidence-stage disease profiles within each RCD concept",
        x=0.02,
        y=0.982,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    fig.text(
        0.02,
        0.065,
        "Radius = evidence-stage signal. Each radar shows the top five disease terms within one death-mode concept after Pair_Count >= 3; disease labels omit death-mode text. "
        "Point size scales with pair PMID count. This is an evidence-stage audit view, not clinical-readiness evidence.",
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#555555",
    )
    fig.subplots_adjust(left=0.035, right=0.80, top=0.89, bottom=0.13, wspace=0.42, hspace=0.52)
    savefig(fig, STEM)

    qc = pd.DataFrame(
        [
            {"metric": "source_table", "value": str(SOURCE)},
            {"metric": "output_source_data", "value": str(out_source)},
            {"metric": "death_modes_plotted", "value": df["death_mode"].nunique()},
            {"metric": "disease_terms_plotted", "value": df["disease_term"].nunique()},
            {"metric": "points_plotted", "value": len(df)},
            {"metric": "radial_metric", "value": "evidence_stage_signal"},
            {"metric": "radial_max", "value": rmax},
            {"metric": "color_source", "value": "Figure 2 death-mode palette"},
            {"metric": "original_figure_preserved", "value": "yes"},
        ]
    )
    qc.to_csv(SRC / f"{STEM}_QC.csv", index=False)
    print(f"Wrote {SUPP / 'pdf' / f'{STEM}.pdf'}")
    print(f"Wrote {SUPP / 'png' / f'{STEM}.png'}")
    print(f"Wrote {SUPP / 'tiff' / f'{STEM}.tiff'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
