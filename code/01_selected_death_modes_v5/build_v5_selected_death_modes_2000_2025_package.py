#!/usr/bin/env python3
"""Build a selected death-mode, 2000-2025 figure package from existing atlas outputs.

This script does not query PubMed and does not change the original v4 package.
It reuses PMID-level and table-level source data from the existing project,
filters all article-level analyses to:

    death_mode in the eight selected regulated cell death modes
    2000 <= publication_year <= 2025

and regenerates the main and supplementary figures in a v4-like package layout.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import textwrap
from datetime import date
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_selected_death_modes")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V4 = Path(os.environ.get("CELLDEATH_ATLAS_V4_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_V4_PACKAGE", "v4_regenerated_source_package")))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
START_YEAR = 2000
END_YEAR = 2025
TODAY = date.today().isoformat()

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

SYSTEM_COLORS = {
    "Neoplasms": "#B9553C",
    "Immune and inflammatory": "#D28A2E",
    "Inflammatory or immune": "#D28A2E",
    "Respiratory": "#3E7CB1",
    "Cardiovascular": "#5C8F7B",
    "Nervous system": "#6B6BAE",
    "Metabolic and nutritional": "#9C6A3D",
    "Wounds and injuries": "#7A9E3A",
    "Musculoskeletal": "#8D6E63",
    "Digestive": "#A05B8F",
    "Digestive system": "#A05B8F",
    "Pathologic process or sign": "#6C9EAE",
    "Other or mixed": "#999999",
}

ONCOLOGY_TERMS = re.compile(
    r"neoplasm|carcinoma|cancer|tumou?r|melanoma|glioma|leukemia|lymphoma|sarcoma|glioblastoma|adenocarcinoma|hepatocellular|myeloma",
    re.I,
)

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


def mkdirs() -> None:
    for d in [
        "00_filter_contract",
        "01_source_data_filtered",
        "02_scripts",
        "03_main_figures_regenerated/pdf",
        "03_main_figures_regenerated/png",
        "03_main_figures_regenerated/tiff",
        "03_main_figures_regenerated/source_data",
        "04_supplementary_figures_regenerated/pdf",
        "04_supplementary_figures_regenerated/png",
        "04_supplementary_figures_regenerated/tiff",
        "04_supplementary_figures_regenerated/source_data",
        "05_manuscript_patches",
        "10_final_qc",
    ]:
        (V5 / d).mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def savefig(fig: plt.Figure, stem: str, main: bool = True) -> None:
    base = V5 / ("03_main_figures_regenerated" if main else "04_supplementary_figures_regenerated")
    for fmt in ["pdf", "png", "tiff"]:
        out = base / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def panel_label(ax, label: str) -> None:
    ax.text(-0.06, 1.06, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def clean_axes(ax, axis: str = "x") -> None:
    ax.grid(axis=axis, color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)


def add_note(ax, text: str, xy=(0.02, 0.02)) -> None:
    ax.text(
        xy[0],
        xy[1],
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#555555",
        bbox=dict(boxstyle="round,pad=0.25", fc="#F7F7F7", ec="#DDDDDD", lw=0.4),
    )


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce").fillna(0.0)
    sd = s.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / sd


def normalize_mode(x: str) -> str:
    if str(x).strip().upper() == "ICD":
        return "Immunogenic cell death"
    return str(x).strip()


def filtered_tmi() -> pd.DataFrame:
    path = V4 / "04_supplementary_figures_regenerated/source_data/tmi_article_stage_audit.csv"
    df = read_csv(path)
    if df.empty:
        raise SystemExit(f"Missing PMID-level audit source: {path}")
    df["death_mode"] = df["death_mode"].map(normalize_mode)
    df["publication_year_num"] = pd.to_numeric(df["publication_year"], errors="coerce")
    df = df[
        df["death_mode"].isin(DEATH_ORDER)
        & df["publication_year_num"].between(START_YEAR, END_YEAR, inclusive="both")
    ].copy()
    df["pmid"] = df["pmid"].astype(str)
    df["disease_term"] = df["disease_term"].astype(str)
    return df


def pair_metrics_base() -> pd.DataFrame:
    # Use the full pair-metrics table, not the v4 Figure 3 display subset.
    # The display subset only contains selected diseases and would misclassify
    # most filtered pairs as "Other or mixed".
    path = Path(
        os.environ.get(
            "CELLDEATH_PAIR_METRICS_BASE",
            PROJECT / "source_data" / "pair_metrics_2025_cutoff.csv",
        )
    )
    df = read_csv(path)
    if df.empty:
        raise SystemExit(f"Missing pair metrics source: {path}")
    df["Death_Mode"] = df["Death_Mode"].map(normalize_mode)
    df = df[df["Death_Mode"].isin(DEATH_ORDER)].copy()
    return df


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").fillna(0.0)
    weights = pd.to_numeric(weights, errors="coerce").fillna(0.0)
    if weights.sum() <= 0:
        return float(values.mean()) if len(values) else 0.0
    return float(np.average(values, weights=weights))


def compute_sources() -> dict[str, pd.DataFrame]:
    tmi = filtered_tmi()
    pair_counts = (
        tmi.groupby(["death_mode", "disease_term"], as_index=False)
        .agg(Pair_Count=("pmid", "nunique"), Pair_PMID_Links=("pmid", "count"))
        .rename(columns={"death_mode": "Death_Mode", "disease_term": "Disease"})
    )

    base = pair_metrics_base()
    base_cols = [
        "Death_Mode",
        "Disease",
        "Disease_Category",
        "Disease_Count",
        "Death_Count",
        "Mechanistic_Convergence_Score",
        "Disease_Specificity_Score",
        "Gap_Score",
        "Translational_Maturity_Index",
        "Article_Count",
    ]
    for col in base_cols:
        if col not in base.columns:
            base[col] = np.nan
    pair = pair_counts.merge(base[base_cols].drop_duplicates(["Death_Mode", "Disease"]), on=["Death_Mode", "Disease"], how="left")
    pair["Disease_Category"] = pair["Disease_Category"].fillna("Other or mixed")
    pair = pair[pair["Pair_Count"] > 0].copy()

    death_total = pair.groupby("Death_Mode")["Pair_Count"].transform("sum")
    disease_total = pair.groupby("Disease")["Pair_Count"].transform("sum")
    total = pair["Pair_Count"].sum()
    pair["Share_Within_Death_Cutoff"] = pair["Pair_Count"] / death_total.replace(0, np.nan)
    pair["Share_Within_Disease_Cutoff"] = pair["Pair_Count"] / disease_total.replace(0, np.nan)
    pair["Disease_Background_Share"] = disease_total / total if total else 0
    pair["Log2_Selectivity_Cutoff"] = np.log2(
        (pair["Share_Within_Death_Cutoff"].fillna(0) + 1e-9)
        / (pair["Disease_Background_Share"].fillna(0) + 1e-9)
    )
    pair["Death_Selectivity_Score"] = pair["Log2_Selectivity_Cutoff"]
    pair["Pair_ID"] = pair["Death_Mode"] + "|||" + pair["Disease"]

    death_summary = (
        pair.groupby("Death_Mode", as_index=False)
        .agg(
            Disease_Breadth=("Disease", "nunique"),
            Total_Pair_Count=("Pair_Count", "sum"),
            Positive_Pairs=("Disease", "size"),
        )
        .merge(
            tmi.groupby("death_mode", as_index=False).agg(Death_Article_N=("pmid", "nunique")).rename(columns={"death_mode": "Death_Mode"}),
            on="Death_Mode",
            how="left",
        )
    )
    death_summary["Death_Mode"] = pd.Categorical(death_summary["Death_Mode"], categories=DEATH_ORDER, ordered=True)
    death_summary = death_summary.sort_values("Death_Mode")

    yearly_rows = []
    for (mode, year), sub in tmi.groupby(["death_mode", "publication_year_num"]):
        dcounts = sub.groupby("disease_term")["pmid"].nunique()
        p = dcounts / dcounts.sum() if dcounts.sum() else dcounts
        entropy = float(-(p * np.log2(p + 1e-12)).sum()) if len(p) else 0.0
        yearly_rows.append(
            {
                "Death_Mode": mode,
                "Year": int(year),
                "Publication_Count": sub["pmid"].nunique(),
                "Active_Disease_Count": sub["disease_term"].nunique(),
                "Disease_Entropy": entropy,
            }
        )
    yearly = pd.DataFrame(yearly_rows)
    all_years = pd.MultiIndex.from_product([DEATH_ORDER, range(START_YEAR, END_YEAR + 1)], names=["Death_Mode", "Year"]).to_frame(index=False)
    yearly = all_years.merge(yearly, on=["Death_Mode", "Year"], how="left").fillna(
        {"Publication_Count": 0, "Active_Disease_Count": 0, "Disease_Entropy": 0}
    )
    summary_rows = []
    for mode in DEATH_ORDER:
        sub = yearly[yearly["Death_Mode"] == mode].copy()
        positive = sub[sub["Publication_Count"] > 0]
        intro = int(positive["Year"].min()) if len(positive) else START_YEAR
        peak = int(sub.loc[sub["Publication_Count"].idxmax(), "Year"]) if len(sub) else END_YEAR
        max_count = sub["Publication_Count"].max()
        take = sub[sub["Publication_Count"] >= max(3, 0.1 * max_count)]
        takeoff = int(take["Year"].min()) if len(take) else intro
        summary_rows.append(
            {
                "Death_Mode": mode,
                "Introduction_Year": intro,
                "Takeoff_Year": takeoff,
                "Peak_Year": peak,
                "Total_Publications": int(positive["Publication_Count"].sum()),
                "Max_Annual_Count": int(max_count),
            }
        )
    temporal_summary = pd.DataFrame(summary_rows)

    hub = compute_hubs(pair)
    oncology = classify_oncology(pair)
    no_neoplasm_pair = pair.merge(oncology[["Disease", "oncology_flag"]], on="Disease", how="left")
    no_neoplasm_pair = no_neoplasm_pair[~no_neoplasm_pair["oncology_flag"].fillna(False)].copy()
    hub_no_neoplasm = compute_hubs(no_neoplasm_pair, score_col="hub_score_no_neoplasms")
    rank_stats = compare_ranks(hub, hub_no_neoplasm)

    mcs_matrix = read_csv(V4 / "03_main_figures_regenerated/source_data/cleaned_death_keyword_log2_rca_matrix.csv")
    mcs_matrix["Death_Mode"] = mcs_matrix["Death_Mode"].map(normalize_mode)
    mcs_matrix = mcs_matrix[mcs_matrix["Death_Mode"].isin(DEATH_ORDER)].copy()

    return {
        "tmi": tmi,
        "pair": pair,
        "death_summary": death_summary,
        "yearly": yearly,
        "temporal_summary": temporal_summary,
        "hub": hub,
        "hub_no_neoplasm": hub_no_neoplasm,
        "oncology": oncology,
        "rank_stats": rank_stats,
        "mcs_matrix": mcs_matrix,
    }


def compute_hubs(pair: pd.DataFrame, score_col: str = "hub_score") -> pd.DataFrame:
    if pair.empty:
        return pd.DataFrame()
    rows = []
    for disease, sub in pair.groupby("Disease"):
        weights = sub["Pair_Count"]
        rows.append(
            {
                "disease_term": disease,
                "disease_system": sub["Disease_Category"].mode().iloc[0] if len(sub["Disease_Category"].dropna()) else "Other or mixed",
                "death_mode_breadth": sub["Death_Mode"].nunique(),
                "pair_count": int(sub["Pair_Count"].sum()),
                "positive_pair_n": len(sub),
                "mean_cleaned_mcs_proxy": weighted_mean(sub["Mechanistic_Convergence_Score"], weights),
                "mean_dss": weighted_mean(sub["Disease_Specificity_Score"], weights),
                "mean_frontier_proxy": weighted_mean(sub["Gap_Score"], weights),
                "unique_pmid_count_proxy": int(sub["Pair_Count"].sum()),
            }
        )
    out = pd.DataFrame(rows)
    out["oncology_flag"] = out.apply(lambda r: is_oncology(r["disease_term"], r["disease_system"]), axis=1)
    out[score_col] = (
        zscore(np.log1p(out["pair_count"]))
        + zscore(out["death_mode_breadth"])
        + zscore(out["mean_cleaned_mcs_proxy"])
        + zscore(out["mean_dss"])
    )
    out = out.sort_values(score_col, ascending=False).reset_index(drop=True)
    out["hub_rank"] = np.arange(1, len(out) + 1)
    if score_col != "hub_score":
        out["hub_score"] = out[score_col]
    return out


def is_oncology(term: str, system: str) -> bool:
    return str(system) == "Neoplasms" or bool(ONCOLOGY_TERMS.search(str(term)))


def classify_oncology(pair: pd.DataFrame) -> pd.DataFrame:
    out = (
        pair.groupby(["Disease", "Disease_Category"], as_index=False)
        .agg(background_count=("Disease_Count", "max"), pair_count_total=("Pair_Count", "sum"))
    )
    out["oncology_flag"] = out.apply(lambda r: is_oncology(r["Disease"], r["Disease_Category"]), axis=1)
    out["oncology_rule_fired"] = np.where(
        out["Disease_Category"].eq("Neoplasms"),
        "disease_system_neoplasms",
        np.where(out["Disease"].str.contains(ONCOLOGY_TERMS, na=False), "term_regex", "not_oncology"),
    )
    return out


def compare_ranks(orig: pd.DataFrame, no: pd.DataFrame) -> pd.DataFrame:
    a = orig[["disease_term", "hub_rank"]].rename(columns={"hub_rank": "rank_original"})
    b = no[["disease_term", "hub_rank"]].rename(columns={"hub_rank": "rank_no_neoplasms"})
    inter = a.merge(b, on="disease_term", how="inner")
    spearman = float(inter["rank_original"].corr(inter["rank_no_neoplasms"], method="spearman")) if len(inter) > 2 else np.nan
    top25_a = set(orig.head(25)["disease_term"])
    top25_b = set(no.head(25)["disease_term"])
    jacc = len(top25_a & top25_b) / len(top25_a | top25_b) if top25_a | top25_b else np.nan
    return pd.DataFrame(
        [
            {
                "comparison": "original_vs_no_neoplasms_selected_death_modes_2000_2025",
                "intersection_entity_n": len(inter),
                "spearman": spearman,
                "top25_jaccard": jacc,
            }
        ]
    )


def save_sources(data: dict[str, pd.DataFrame]) -> None:
    mapping = {
        "tmi": "article_stage_records_selected_death_modes_2000_2025.csv",
        "pair": "pair_metrics_selected_death_modes_2000_2025.csv",
        "death_summary": "death_summary_selected_death_modes_2000_2025.csv",
        "yearly": "death_yearly_selected_death_modes_2000_2025.csv",
        "temporal_summary": "temporal_summary_selected_death_modes_2000_2025.csv",
        "hub": "hub_score_original_selected_death_modes_2000_2025.csv",
        "hub_no_neoplasm": "hub_score_no_neoplasms_selected_death_modes_2000_2025.csv",
        "oncology": "oncology_term_classification_selected_death_modes_2000_2025.csv",
        "rank_stats": "oncology_rank_stats_selected_death_modes_2000_2025.csv",
        "mcs_matrix": "cleaned_death_keyword_log2_rca_matrix_selected_death_modes.csv",
    }
    for key, name in mapping.items():
        df = data[key]
        df.to_csv(V5 / "01_source_data_filtered" / name, index=False)
        target = V5 / ("03_main_figures_regenerated/source_data" if key in {"pair", "death_summary", "yearly", "temporal_summary", "hub", "hub_no_neoplasm", "mcs_matrix", "rank_stats"} else "04_supplementary_figures_regenerated/source_data") / name
        target.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(target, index=False)


def build_filter_contract(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"]
    tmi = data["tmi"]
    text = f"""
    # Eight-Concept 2000-2025 Figure Regeneration Contract

    ## Scope

    - Included concepts: {', '.join(DEATH_ORDER)}
    - Publication-year window: {START_YEAR}-{END_YEAR}
    - PubMed is not requeried.
    - Existing v4 remains unchanged.

    ## Regenerated atlas scale after filtering

    - Death-mode concepts: {len(DEATH_ORDER)}
    - Disease-system mapped terms with at least one filtered pair: {pair['Disease'].nunique():,}
    - Positive death-mode disease pairs: {len(pair):,}
    - Pair-PMID article-stage rows: {len(tmi):,}
    - Unique PubMed records: {tmi['pmid'].nunique():,}

    ## Interpretation

    The regenerated figures describe an eight-concept 2000-2025 literature atlas. All counts differ from v4 because the final selected death-mode scope and publication-year window are applied. Co-mention remains literature association, not causal evidence.
    """
    write_text(V5 / "00_filter_contract/selected_death_modes_2000_2025_contract.md", text)


def draw_workflow(ax) -> None:
    ax.axis("off")
    nodes = [
        ("Death-mode query definitions\n8 ontology-qualified concepts", 0.02, 0.66, "#EEF4F8"),
        ("MeSH-derived disease-term set\ncleaned and system mapped", 0.02, 0.30, "#F4F1E8"),
        ("Pairwise PubMed retrieval\ndeath-mode x disease queries", 0.33, 0.48, "#F7F7F7"),
        ("Filter retained records\n2000-2025 retained records", 0.54, 0.48, "#FFF4E6"),
        ("Pair-PMID link table\narticle-level associations", 0.72, 0.48, "#F7F7F7"),
        ("Audited atlas layers\nselectivity, hubs, cleaned vocabulary,\nevidence-stage audit, robustness", 0.72, 0.08, "#EAF3EF"),
    ]
    for text, x, y, color in nodes:
        w = 0.18 if x < 0.54 else 0.16
        box = FancyBboxPatch((x, y), w, 0.18, boxstyle="round,pad=0.012,rounding_size=0.02", fc=color, ec="#444444", lw=0.65)
        ax.add_patch(box)
        ax.text(x + w / 2, y + 0.09, text, ha="center", va="center", fontsize=6.0)
    arrows = [((0.20, 0.75), (0.33, 0.57)), ((0.20, 0.39), (0.33, 0.52)), ((0.51, 0.57), (0.54, 0.57)), ((0.70, 0.57), (0.72, 0.57)), ((0.80, 0.48), (0.80, 0.26))]
    for start, end in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#555555"))
    ax.text(0.02, 0.02, "Co-mention indicates literature association, not causal evidence.", fontsize=6.2, color="#555555")


def figure1(data: dict[str, pd.DataFrame]) -> None:
    death = data["death_summary"].copy()
    pair = data["pair"].copy()
    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1], width_ratios=[1.1, 0.9])
    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    panel_label(ax0, "A")
    draw_workflow(ax0)
    panel_label(ax1, "B")
    death = death.set_index("Death_Mode").reindex(DEATH_ORDER).reset_index()
    ax1.scatter(
        death["Death_Article_N"],
        death["Disease_Breadth"],
        s=np.clip(np.sqrt(death["Total_Pair_Count"]) * 10, 60, 620),
        c=[MODE_COLORS[m] for m in death["Death_Mode"]],
        alpha=0.86,
        edgecolor="white",
        linewidth=0.7,
    )
    label_offsets = {
        "Ferroptosis": (10, 8),
        "Pyroptosis": (8, -11),
        "NETosis": (-38, 9),
        "Necroptosis": (8, -2),
        "Immunogenic cell death": (8, -10),
        "Cuproptosis": (8, 8),
        "PANoptosis": (8, -10),
        "Disulfidptosis": (8, -2),
    }
    for _, r in death.iterrows():
        dx, dy = label_offsets.get(r["Death_Mode"], (8, 0))
        ha = "right" if dx < 0 else "left"
        ax1.annotate(
            r["Death_Mode"],
            (r["Death_Article_N"], r["Disease_Breadth"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.5,
            ha=ha,
            va="center",
        )
    ax1.set_xscale("log")
    ax1.set_xlabel("Death-mode unique PMIDs, log scale")
    ax1.set_ylabel("Disease breadth")
    ax1.set_title("Uneven breadth and volume after 2000-2025 filtering", loc="left", weight="bold")
    clean_axes(ax1)
    panel_label(ax2, "C")
    topcat = pair["Disease_Category"].value_counts().head(8)
    ax2.barh(range(len(topcat))[::-1], topcat.values, color=[SYSTEM_COLORS.get(k, "#BDBDBD") for k in topcat.index], height=0.72)
    ax2.set_yticks(range(len(topcat))[::-1])
    ax2.set_yticklabels(topcat.index)
    ax2.set_xlabel("Positive pairs")
    ax2.set_title("Disease-system composition", loc="left", weight="bold")
    clean_axes(ax2)
    fig.suptitle("Figure 1. Eight-concept PubMed atlas construction and disease-literature scale", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Figure_1_final_selected_death_modes_2000_2025", True)


def figure2(data: dict[str, pd.DataFrame]) -> None:
    yearly = data["yearly"]
    summary = data["temporal_summary"]
    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.25, 0.75])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    panel_label(ax0, "A")
    eras = [(2000, 2006, "early disease-facing\nRCD vocabulary"), (2007, 2014, "inflammatory death\nadoption"), (2015, 2020, "ferroptosis/redox\nexpansion"), (2021, 2025, "metal-stress and\nintegrative phase")]
    era_cols = ["#F4F1E8", "#F8EBDD", "#EEF4F8", "#EAF3EF"]
    for (x0, x1, label), col in zip(eras, era_cols):
        ax0.axvspan(x0, x1, color=col, zorder=0)
        ax0.text((x0 + x1) / 2, 1.08, label, ha="center", va="bottom", fontsize=6.0, color="#555555", transform=ax0.get_xaxis_transform())
    for mode in DEATH_ORDER:
        sub = yearly[yearly["Death_Mode"] == mode].sort_values("Year")
        y = sub["Publication_Count"].astype(float)
        ax0.plot(sub["Year"], y, lw=1.45, color=MODE_COLORS[mode], label=mode)
    ax0.set_xlim(START_YEAR, END_YEAR)
    max_y = float(pd.to_numeric(yearly["Publication_Count"], errors="coerce").max())
    ax0.set_ylim(0, max_y * 1.16 if max_y > 0 else 1)
    ax0.set_ylabel("Annual publication count, unique PMIDs")
    ax0.set_title("Annual publication volume in disease literature", loc="left", weight="bold")
    ax0.legend(ncol=3, bbox_to_anchor=(0, -0.05), loc="upper left", frameon=False)
    clean_axes(ax0)
    panel_label(ax1, "B")
    summary = summary.set_index("Death_Mode").reindex(DEATH_ORDER).reset_index()
    y = np.arange(len(summary))[::-1]
    ax1.hlines(y, summary["Introduction_Year"], summary["Peak_Year"], color="#D0D0D0", lw=2)
    ax1.scatter(summary["Introduction_Year"], y, s=28, color="#FFFFFF", edgecolor="#555555", zorder=3, label="first year in filtered atlas")
    ax1.scatter(summary["Takeoff_Year"], y, s=38, color=[MODE_COLORS[m] for m in summary["Death_Mode"]], zorder=3, label="takeoff")
    ax1.scatter(summary["Peak_Year"], y, s=24, marker="s", color="#333333", zorder=3, label="peak")
    ax1.set_yticks(y)
    ax1.set_yticklabels(summary["Death_Mode"])
    ax1.set_xlim(START_YEAR, END_YEAR + 1)
    ax1.set_xlabel("Year")
    ax1.set_title("Takeoff years after 2000-2025 filtering", loc="left", weight="bold")
    ax1.legend(ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, borderaxespad=0)
    clean_axes(ax1)
    fig.suptitle("Figure 2. Historical waves of RCD terminology adoption, 2000-2025", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Figure_2_final_selected_death_modes_2000_2025", True)


def figure3(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"].copy()
    focus = [
        ("Neoplasms", 5),
        ("Immune and inflammatory", 3),
        ("Respiratory", 3),
        ("Cardiovascular", 3),
        ("Metabolic and nutritional", 3),
        ("Nervous system", 3),
        ("Digestive", 3),
        ("Musculoskeletal", 3),
        ("Pathologic process or sign", 3),
    ]
    diseases: list[str] = []
    for cat, n in focus:
        sub = pair[pair["Disease_Category"].eq(cat)].groupby(["Disease", "Disease_Category"], as_index=False)["Pair_Count"].sum().sort_values("Pair_Count", ascending=False).head(n)
        diseases.extend(sub["Disease"].tolist())
    diseases = list(dict.fromkeys(diseases))
    d = pair[pair["Disease"].isin(diseases)].copy()
    order = d.groupby(["Disease_Category", "Disease"])["Pair_Count"].sum().reset_index().sort_values(["Disease_Category", "Pair_Count"], ascending=[True, False])
    diseases = order["Disease"].drop_duplicates().tolist()[::-1]
    xmap = {m: i for i, m in enumerate(DEATH_ORDER)}
    ymap = {disease: i for i, disease in enumerate(diseases)}
    d["x"] = d["Death_Mode"].map(xmap)
    d["y"] = d["Disease"].map(ymap)
    fig, ax = plt.subplots(figsize=(7.4, max(5.8, 0.22 * len(diseases))), constrained_layout=True)
    panel_label(ax, "A")
    sizes = np.clip(np.sqrt(d["Pair_Count"].astype(float)) * 9, 12, 240)
    sc = ax.scatter(d["x"], d["y"], s=sizes, c=d["Log2_Selectivity_Cutoff"], cmap="viridis", alpha=0.84, edgecolor="white", lw=0.35)
    ax.set_xticks(range(len(DEATH_ORDER)))
    ax.set_xticklabels(DEATH_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(diseases)))
    ax.set_yticklabels(diseases)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.62, pad=0.01)
    cbar.set_label("Disease-literature selectivity")
    for cat, sub in order.groupby("Disease_Category", sort=False):
        ys = [ymap[x] for x in sub["Disease"] if x in ymap]
        if ys:
            ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color=SYSTEM_COLORS.get(cat, "#EEEEEE"), alpha=0.08, zorder=-1)
            ax.text(len(DEATH_ORDER) + 0.12, np.mean(ys), cat, va="center", fontsize=5.8, color=SYSTEM_COLORS.get(cat, "#555555"))
    ax.set_xlim(-0.6, len(DEATH_ORDER) + 2.8)
    ax.grid(color="#EAEAEA", lw=0.35)
    ax.set_title("Disease niches after 2000-2025 filtering", loc="left", weight="bold")
    add_note(ax, "Bubble area: 2000-2025 pair count; color: disease-literature selectivity.")
    fig.suptitle("Figure 3. Disease-literature selectivity across eight RCD-associated concepts", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Figure_3_final_selected_death_modes_2000_2025", True)


def figure4(data: dict[str, pd.DataFrame]) -> None:
    orig = data["hub"].head(16).copy()
    no = data["hub_no_neoplasm"].head(12).copy()
    stats = data["rank_stats"].iloc[0].to_dict()
    oncology = data["oncology"]
    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.28])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])
    panel_label(ax0, "A")
    y0 = np.arange(len(orig))[::-1]
    ax0.barh(y0, orig["hub_score"], color=["#B9553C" if x else "#6C9EAE" for x in orig["oncology_flag"]], height=0.68)
    ax0.set_yticks(y0)
    ax0.set_yticklabels(orig["disease_term"])
    ax0.set_xlabel("Hub score")
    ax0.set_title("Original eight-concept hubs", loc="left", weight="bold")
    clean_axes(ax0)
    panel_label(ax1, "B")
    y1 = np.arange(len(no))[::-1]
    ax1.barh(y1, no["hub_score_no_neoplasms"], color=[SYSTEM_COLORS.get(x, "#6C9EAE") for x in no["disease_system"]], height=0.68)
    ax1.set_yticks(y1)
    ax1.set_yticklabels(no["disease_term"])
    ax1.set_xlabel("Hub score after neoplasm removal")
    ax1.set_title("Retained non-oncology structure", loc="left", weight="bold")
    clean_axes(ax1)
    panel_label(ax2, "C")
    ax2.axis("off")
    cards = [
        (f"{int(oncology['oncology_flag'].sum()):,}", "oncology/neoplasm\nterms flagged", "#F7EAE6"),
        (f"{int((~oncology['oncology_flag']).sum()):,}", "non-oncology\nterms retained", "#EAF3EF"),
        (f"Spearman = {stats['spearman']:.3f}", "global ranking\nstructure", "#EEF4F8"),
        (f"Top25 Jaccard = {stats['top25_jaccard']:.3f}", "top hub identity\noverlap", "#F4F1E8"),
    ]
    for i, (big, small, color) in enumerate(cards):
        x = 0.02 + i * 0.245
        ax2.add_patch(FancyBboxPatch((x, 0.20), 0.22, 0.58, boxstyle="round,pad=0.012,rounding_size=0.02", fc=color, ec="#CCCCCC", lw=0.6))
        ax2.text(x + 0.11, 0.57, big, ha="center", va="center", weight="bold", fontsize=9)
        ax2.text(x + 0.11, 0.35, small, ha="center", va="center", fontsize=6.3)
    fig.suptitle("Figure 4. Oncology-sensitive literature hubs in the eight-concept atlas", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Figure_4_final_selected_death_modes_2000_2025", True)


def figure5(data: dict[str, pd.DataFrame]) -> None:
    mat = data["mcs_matrix"].set_index("Death_Mode")
    axis_terms = {
        "metabolic-redox": ["gpx4", "acsl4", "lipid peroxidation", "iron overload", "ferritinophagy", "ncoa4", "mitochondrial dysfunction", "oxidative stress", "mitophagy", "metabolic reprogramming"],
        "inflammatory/innate": ["inflammasome", "gasdermin d", "gasdermin", "il-1b", "il-6", "nf-kb", "interferon", "chemokine", "macrophage polarization"],
        "tumor-immune": ["pd-l1", "t cell exhaustion", "pd-1", "damage-associated molecular pattern", "ifn-gamma"],
        "integrative death": ["ripk1", "ripk3", "mlkl", "caspase 8", "fadd", "caspase 1", "caspase 3"],
    }
    selected, groups = [], []
    cols_lower = {c.lower(): c for c in mat.columns}
    forbidden = {"machine learning", "mice, knockout", "knockout mice"}
    for group, terms in axis_terms.items():
        for term in terms:
            if term in cols_lower and term not in forbidden and cols_lower[term] not in selected:
                selected.append(cols_lower[term])
                groups.append(group)
    data_mat = mat.reindex(DEATH_ORDER)[selected].fillna(0)
    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 0.12])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    panel_label(ax, "A")
    sns.heatmap(data_mat, ax=ax, cmap="vlag", center=0, linewidths=0.25, linecolor="#FFFFFF", cbar_kws={"label": "log2 RCA", "shrink": 0.55})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Cleaned mechanism-vocabulary axes in the eight-concept atlas", loc="left", weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax2.axis("off")
    palette = {"metabolic-redox": "#B9553C", "inflammatory/innate": "#D28A2E", "tumor-immune": "#A05B8F", "integrative death": "#5C8F7B"}
    x = 0
    for group in ["metabolic-redox", "inflammatory/innate", "tumor-immune", "integrative death"]:
        n = groups.count(group)
        if n:
            ax2.add_patch(Rectangle((x / len(selected), 0.2), n / len(selected), 0.45, transform=ax2.transAxes, color=palette[group], alpha=0.75))
            ax2.text((x + n / 2) / len(selected), 0.78, group, transform=ax2.transAxes, ha="center", va="bottom", fontsize=6.1)
            x += n
    add_note(ax, "Vocabulary-level signals; not experimental proof of pathway equivalence.")
    fig.suptitle("Figure 5. Cleaned mechanism-vocabulary axes across eight-concept RCD-associated literature", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Figure_5_final_selected_death_modes_2000_2025", True)


def threshold_summary(pair: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for thr in [1, 3, 5, 10]:
        sub = pair[pair["Pair_Count"] >= thr]
        rows.append(
            {
                "threshold": thr,
                "threshold_label": f"Pair_Count >= {thr}",
                "retained_pair_count": len(sub),
                "retained_disease_count": sub["Disease"].nunique(),
                "retained_death_mode_breadth": sub["Death_Mode"].nunique(),
                "retained_pair_count_volume": int(sub["Pair_Count"].sum()),
            }
        )
    return pd.DataFrame(rows)


def s1(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"]
    thr = threshold_summary(pair)
    thr.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/threshold_summary_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    panel_label(axes[0], "A")
    axes[0].hist(np.clip(pair["Pair_Count"], 1, 50), bins=np.arange(1, 52), color="#8FB3C8", edgecolor="white")
    for x in [3, 10]:
        axes[0].axvline(x, color="#444444", ls="--", lw=0.8)
        axes[0].text(x + 0.25, axes[0].get_ylim()[1] * 0.82, f">={x}", fontsize=6)
    axes[0].set_xlabel("Pair count, clipped at 50")
    axes[0].set_ylabel("Pair number")
    axes[0].set_title("Low-count relationships remain numerous", loc="left", weight="bold")
    panel_label(axes[1], "B")
    axes[1].plot(thr["threshold"], thr["retained_pair_count"], marker="o", color="#3E7CB1", label="pairs")
    axes[1].plot(thr["threshold"], thr["retained_disease_count"], marker="o", color="#B9553C", label="diseases")
    axes[1].set_xlabel("Pair_Count threshold")
    axes[1].set_ylabel("Retained count")
    axes[1].set_title("Threshold retention", loc="left", weight="bold")
    axes[1].legend(frameon=False)
    for ax in axes:
        clean_axes(ax)
    savefig(fig, "Supplementary_Figure_S1_low_count_structure_selected_death_modes_2000_2025", False)


def rank_stability_tables(pair: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    hubs = {}
    metrics = {}
    for thr in [1, 3, 5, 10]:
        sub = pair[pair["Pair_Count"] >= thr]
        h = compute_hubs(sub)
        hubs[thr] = h[["disease_term", "hub_rank", "hub_score"]].copy()
        m = sub.copy()
        for metric, col in [("cleaned_mcs", "Mechanistic_Convergence_Score"), ("dss", "Disease_Specificity_Score"), ("frontier", "Gap_Score")]:
            rank = (
                m.groupby("Disease", as_index=False)[col]
                .mean()
                .sort_values(col, ascending=False)
                .reset_index(drop=True)
                .rename(columns={"Disease": "entity", col: "value"})
            )
            rank["rank"] = np.arange(1, len(rank) + 1)
            metrics[(metric, thr)] = rank
    rows = []
    metric_rows = []
    ref = hubs[1]
    for thr, h in hubs.items():
        inter = ref.merge(h, on="disease_term", suffixes=("_ref", "_cmp"))
        rows.append(
            {
                "metric_name": "hub_score",
                "reference_threshold": 1,
                "comparison_threshold": thr,
                "intersection_entity_n": len(inter),
                "spearman_correlation": inter["hub_rank_ref"].corr(inter["hub_rank_cmp"], method="spearman") if len(inter) > 2 else np.nan,
                "top_10_jaccard": jaccard(ref.head(10)["disease_term"], h.head(10)["disease_term"]),
                "top_25_jaccard": jaccard(ref.head(25)["disease_term"], h.head(25)["disease_term"]),
                "top_50_jaccard": jaccard(ref.head(50)["disease_term"], h.head(50)["disease_term"]),
            }
        )
    for metric in ["cleaned_mcs", "dss", "frontier"]:
        refm = metrics[(metric, 1)]
        for thr in [1, 3, 5, 10]:
            cmp = metrics[(metric, thr)]
            inter = refm.merge(cmp, on="entity", suffixes=("_ref", "_cmp"))
            metric_rows.append(
                {
                    "metric_name": metric,
                    "reference_threshold": 1,
                    "comparison_threshold": thr,
                    "intersection_entity_n": len(inter),
                    "spearman_correlation": inter["rank_ref"].corr(inter["rank_cmp"], method="spearman") if len(inter) > 2 else np.nan,
                    "top_25_jaccard": jaccard(refm.head(25)["entity"], cmp.head(25)["entity"]),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(metric_rows)


def jaccard(a: Iterable, b: Iterable) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if sa | sb else np.nan


def s2(data: dict[str, pd.DataFrame]) -> None:
    hub_stab, metric_stab = rank_stability_tables(data["pair"])
    hub_stab.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/hub_rank_stability_selected_death_modes_2000_2025.csv", index=False)
    metric_stab.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/metric_rank_stability_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    panel_label(axes[0], "A")
    pivot = hub_stab.pivot_table(index="comparison_threshold", values=["spearman_correlation", "top_25_jaccard"], aggfunc="mean")
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1, ax=axes[0], cbar=False)
    axes[0].set_title("Hub stability across thresholds", loc="left", weight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Threshold")
    panel_label(axes[1], "B")
    m = metric_stab.groupby("metric_name")["spearman_correlation"].mean().sort_values()
    axes[1].barh(m.index, m.values, color="#6C9EAE")
    axes[1].set_xlim(0, 1.02)
    axes[1].set_xlabel("Mean Spearman")
    axes[1].set_title("Metric rank stability", loc="left", weight="bold")
    clean_axes(axes[1])
    savefig(fig, "Supplementary_Figure_S2_rank_stability_selected_death_modes_2000_2025", False)


def s3(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"].copy()
    rng = np.random.default_rng(20260513)
    observed = {
        "mean_pair_count": pair["Pair_Count"].mean(),
        "max_hub_score": data["hub"]["hub_score"].max(),
        "disease_system_pair_variance": pair["Disease_Category"].value_counts(normalize=True).var(),
        "max_frontier_proxy": pair["Gap_Score"].max(skipna=True),
    }
    rows = []
    for metric_name, obs in observed.items():
        vals = []
        for _ in range(200):
            shuf = pair.copy()
            shuf["Disease"] = rng.permutation(shuf["Disease"].values)
            shuf["Disease_Category"] = rng.permutation(shuf["Disease_Category"].values)
            if metric_name == "mean_pair_count":
                vals.append(shuf["Pair_Count"].mean())
            elif metric_name == "max_hub_score":
                vals.append(compute_hubs(shuf)["hub_score"].max())
            elif metric_name == "disease_system_pair_variance":
                vals.append(shuf["Disease_Category"].value_counts(normalize=True).var())
            elif metric_name == "max_frontier_proxy":
                vals.append(shuf["Gap_Score"].max(skipna=True))
        vals = np.array(vals, dtype=float)
        rows.append(
            {
                "null_model": "shuffle_disease_labels",
                "metric_name": metric_name,
                "observed_value": obs,
                "null_mean": np.nanmean(vals),
                "null_sd": np.nanstd(vals),
                "null_ci_lower": np.nanpercentile(vals, 2.5),
                "null_ci_upper": np.nanpercentile(vals, 97.5),
                "permutation_iterations": 200,
                "interpretation_note": "Exploratory null control only; do not claim formal statistical significance.",
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/permutation_null_summary_selected_death_modes_2000_2025.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 3.5), constrained_layout=True)
    panel_label(ax, "A")
    d = df.iloc[::-1].copy()
    y = np.arange(len(d))
    ax.hlines(y, d["null_ci_lower"], d["null_ci_upper"], color="#BDBDBD", lw=2)
    ax.scatter(d["null_mean"], y, color="#888888", s=24, label="null mean")
    ax.scatter(d["observed_value"], y, color="#B9553C", s=28, label="observed")
    ax.set_yticks(y)
    ax.set_yticklabels(d["metric_name"].str.replace("_", " "))
    ax.set_xlabel("Metric value")
    ax.set_title("Observed structure versus shuffled disease-label nulls", loc="left", weight="bold")
    ax.legend(frameon=False, loc="lower right")
    clean_axes(ax)
    add_note(ax, "Exploratory structural diagnostic only; no formal significance claimed.")
    savefig(fig, "Supplementary_Figure_S3_observed_vs_null_selected_death_modes_2000_2025", False)


def s4(data: dict[str, pd.DataFrame]) -> bool:
    for subdir in ["pdf", "png", "tiff"]:
        for old in (V5 / "04_supplementary_figures_regenerated" / subdir).glob("Supplementary_Figure_S4*.pdf" if subdir == "pdf" else "Supplementary_Figure_S4*.*"):
            old.unlink()

    metrics_path = PROJECT / "reviewer_revision/11_supplementary_recovery/journal_atlas_share_metrics_2000_2025_selected_death_modes.csv"
    missing_path = V5 / "04_supplementary_figures_regenerated/Supplementary_Figure_S4_missing_source_report.md"
    if not metrics_path.exists():
        write_text(
            missing_path,
            f"""
            # Supplementary Figure S4 Missing Source

            Journal denominator counts were not found for the strict selected death-mode {START_YEAR}-{END_YEAR} window.

            Required input:

            `{metrics_path}`

            Regenerate this file with:

            `python3 scripts/reviewer_revision/journal_denominator_pubmed_2000_2025_selected_death_modes.py --run-all --email YOUR_EMAIL@example.com`

            The older `journal_atlas_share_metrics.csv` file is intentionally not accepted because it uses a different denominator time window.
            """,
        )
        return False
    met = read_csv(metrics_path)
    required = {
        "journal_display_name",
        "journal_normalized",
        "atlas_unique_pmids",
        "denominator_total_pubmed_articles",
        "atlas_share_percent",
        "death_mode_diversity",
        "denominator_year_min",
        "denominator_year_max",
        "included_in_main_plot_yes_no",
    }
    missing = sorted(required.difference(met.columns))
    if missing:
        write_text(
            missing_path,
            f"""
            # Supplementary Figure S4 Missing Source

            The strict selected death-mode journal metrics file exists but is missing required columns.

            Required input: `{metrics_path}`

            Missing columns: {', '.join(missing)}
            """,
        )
        return False
    met["denominator_year_min"] = pd.to_numeric(met["denominator_year_min"], errors="coerce")
    met["denominator_year_max"] = pd.to_numeric(met["denominator_year_max"], errors="coerce")
    met = met[(met["denominator_year_min"] == START_YEAR) & (met["denominator_year_max"] == END_YEAR)].copy()
    if met.empty:
        write_text(
            missing_path,
            f"""
            # Supplementary Figure S4 Missing Source

            The journal metrics file does not contain rows with denominator_year_min={START_YEAR} and denominator_year_max={END_YEAR}.

            Required input: `{metrics_path}`

            S4 was not regenerated to avoid mixing the selected death-mode {START_YEAR}-{END_YEAR} numerator with a mismatched journal denominator.
            """,
        )
        return False
    met["denominator_total_pubmed_articles"] = pd.to_numeric(met["denominator_total_pubmed_articles"], errors="coerce")
    met["atlas_unique_pmids"] = pd.to_numeric(met["atlas_unique_pmids"], errors="coerce")
    met["atlas_share_percent"] = pd.to_numeric(met["atlas_share_percent"], errors="coerce")
    met["death_mode_diversity"] = pd.to_numeric(met["death_mode_diversity"], errors="coerce")
    met["included_in_plot"] = met["included_in_main_plot_yes_no"].astype(str).str.lower().isin(["yes", "true", "1"])
    met = met.sort_values("atlas_unique_pmids", ascending=False)
    met.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/journal_atlas_share_metrics_selected_death_modes_2000_2025.csv", index=False)
    plot = met[met["included_in_plot"]].copy()
    if plot.empty:
        write_text(missing_path, "# Supplementary Figure S4 Missing Source\n\nNo journals passed denominator and numerator filters after selected death-mode 2000-2025 filtering; S4 was not regenerated.\n")
        return False
    fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])
    panel_label(ax0, "A")
    top_count = plot.head(12).iloc[::-1]
    ax0.barh(top_count["journal_display_name"], top_count["atlas_unique_pmids"], color="#6C9EAE")
    ax0.set_xlabel("Filtered atlas PMIDs")
    ax0.set_title("Top journals by atlas count", loc="left", weight="bold")
    clean_axes(ax0)
    panel_label(ax1, "B")
    top_share = plot.sort_values("atlas_share_percent", ascending=False).head(12).iloc[::-1]
    ax1.barh(top_share["journal_display_name"], top_share["atlas_share_percent"], color="#D28A2E")
    ax1.set_xlabel("Atlas share of PubMed output (%)")
    ax1.set_title("Denominator-normalized journal context", loc="left", weight="bold")
    clean_axes(ax1)
    panel_label(ax2, "C")
    sc = ax2.scatter(
        plot["denominator_total_pubmed_articles"],
        plot["atlas_unique_pmids"],
        s=np.clip(plot["atlas_share_percent"] * 20, 20, 250),
        c=plot["death_mode_diversity"],
        cmap="viridis",
        alpha=0.75,
        edgecolor="white",
    )
    ax2.set_xscale("log")
    ax2.set_xlabel("Total PubMed-indexed journal articles, log scale")
    ax2.set_ylabel("Filtered atlas PMIDs")
    ax2.set_title("Atlas count versus journal denominator", loc="left", weight="bold")
    cbar = fig.colorbar(sc, ax=ax2, shrink=0.6, pad=0.01)
    cbar.set_label("Death-mode diversity")
    clean_axes(ax2)
    add_note(ax2, f"{START_YEAR}-{END_YEAR} denominator-normalized journal context only; not a journal quality or biological-importance metric.")
    fig.suptitle("Supplementary Figure S4. Journal-level context after denominator normalization", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Supplementary_Figure_S4_journal_context_selected_death_modes_2000_2025", False)
    return True


def s5(data: dict[str, pd.DataFrame]) -> None:
    no = data["hub_no_neoplasm"].head(18)
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    panel_label(ax, "A")
    y = np.arange(len(no))[::-1]
    ax.barh(y, no["hub_score_no_neoplasms"], color=[SYSTEM_COLORS.get(x, "#6C9EAE") for x in no["disease_system"]], height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(no["disease_term"])
    ax.set_xlabel("Hub score after neoplasm removal")
    ax.set_title("Retained non-oncology literature hubs after filtering", loc="left", weight="bold")
    clean_axes(ax)
    savefig(fig, "Supplementary_Figure_S5_retained_non_oncology_hubs_selected_death_modes_2000_2025", False)


def s6(data: dict[str, pd.DataFrame]) -> None:
    tmi = data["tmi"]
    valid_guideline = tmi["valid_guideline_yes_no"].astype(str).str.lower().eq("yes").sum() if "valid_guideline_yes_no" in tmi else 0
    guideline_ineligible = tmi["guideline_evidence_status"].astype(str).eq("guideline_ineligible_in_article_corpus").sum() if "guideline_evidence_status" in tmi else 0
    original_guideline = tmi["original_stage"].astype(str).eq("guideline_evidence").sum() if "original_stage" in tmi else guideline_ineligible
    trial_supported = tmi["interventional_supported_by_publication_type_yes_no"].astype(str).str.lower().eq("yes").sum() if "interventional_supported_by_publication_type_yes_no" in tmi else 0
    by_mode = (
        tmi.groupby("death_mode", as_index=False)
        .agg(
            original_high=("original_high_stage_yes_no", lambda x: x.astype(str).str.lower().eq("yes").sum()),
            audited_valid_high=("audited_stage", lambda x: x.astype(str).isin(["interventional_trial", "guideline_evidence"]).sum()),
        )
    )
    by_mode["death_mode"] = pd.Categorical(by_mode["death_mode"], categories=DEATH_ORDER, ordered=True)
    by_mode = by_mode.sort_values("death_mode")
    by_mode.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/tmi_summary_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7), constrained_layout=True)
    panel_label(axes[0], "A")
    vals = [original_guideline, valid_guideline, trial_supported]
    labs = ["guideline-like\noriginal", "valid guideline\nevidence", "PublicationType-\nsupported trials"]
    axes[0].bar(labs, vals, color=["#C7C7C7", "#6C9EAE", "#3E7CB1"])
    axes[0].set_ylabel("Article assignments")
    axes[0].set_title("Guideline signal remains ineligible", loc="left", weight="bold")
    panel_label(axes[1], "B")
    x = np.arange(len(by_mode))
    axes[1].bar(x - 0.18, by_mode["original_high"], width=0.35, color="#C7C7C7", label="original high-stage")
    axes[1].bar(x + 0.18, by_mode["audited_valid_high"], width=0.35, color="#3E7CB1", label="audited valid high-stage")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(by_mode["death_mode"], rotation=45, ha="right")
    axes[1].set_ylabel("Assignments")
    axes[1].set_title("Evidence-stage audit by concept", loc="left", weight="bold")
    axes[1].legend(frameon=False)
    for ax in axes:
        clean_axes(ax)
    savefig(fig, "Supplementary_Figure_S6_evidence_stage_audit_selected_death_modes_2000_2025", False)


def s7(data: dict[str, pd.DataFrame]) -> None:
    samp_path = V4 / "04_supplementary_figures_regenerated/source_data/manual_curation_sample.csv"
    samp = read_csv(samp_path)
    if not samp.empty and "death_mode" in samp:
        samp = samp[samp["death_mode"].map(normalize_mode).isin(DEATH_ORDER)].copy()
    samp.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/manual_curation_sample_selected_death_modes_subset.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    panel_label(axes[0], "A")
    axes[0].axis("off")
    steps = ["Filtered\nPMID sample", "Dual human\ncuration", "Labels\nblank", "Agreement\nfuture step"]
    for i, step in enumerate(steps):
        x = 0.03 + i * 0.24
        axes[0].add_patch(FancyBboxPatch((x, 0.43), 0.18, 0.22, boxstyle="round,pad=0.02", fc="#F7F7F7", ec="#777777", lw=0.6))
        axes[0].text(x + 0.09, 0.54, step, ha="center", va="center", fontsize=6.5)
        if i < len(steps) - 1:
            axes[0].add_patch(FancyArrowPatch((x + 0.19, 0.54), (x + 0.235, 0.54), arrowstyle="-|>", mutation_scale=8, lw=0.7, color="#555555"))
    axes[0].set_title("Manual-curation framework, not completed validation", loc="left", weight="bold")
    panel_label(axes[1], "B")
    if not samp.empty and "stratum" in samp:
        counts = samp["stratum"].value_counts().head(12).iloc[::-1]
        axes[1].barh(counts.index, counts.values, color="#8FB3C8")
    axes[1].set_xlabel("Sampled PMIDs after filtering")
    axes[1].set_title("Remaining curation strata", loc="left", weight="bold")
    clean_axes(axes[1])
    savefig(fig, "Supplementary_Figure_S7_manual_curation_framework_selected_death_modes_2000_2025", False)


def s8(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"]
    # High-risk low-count pairs are recomputed conservatively from filtered data.
    mcs_cut = pair["Mechanistic_Convergence_Score"].quantile(0.9)
    tmi_cut = pair["Translational_Maturity_Index"].quantile(0.9)
    frontier_cut = pair["Gap_Score"].quantile(0.9)
    rows = []
    for _, r in pair[pair["Pair_Count"] < 5].iterrows():
        reasons = []
        if pd.notna(r["Mechanistic_Convergence_Score"]) and r["Mechanistic_Convergence_Score"] >= mcs_cut:
            reasons.append("pair_count_lt5_and_high_cleaned_MCS")
        if pd.notna(r["Translational_Maturity_Index"]) and r["Translational_Maturity_Index"] >= tmi_cut:
            reasons.append("pair_count_lt5_and_high_clinical_signal")
        if pd.notna(r["Gap_Score"]) and r["Gap_Score"] >= frontier_cut:
            reasons.append("pair_count_lt5_and_high_Frontier")
        if reasons:
            rows.append(
                {
                    "death_mode": r["Death_Mode"],
                    "disease_term": r["Disease"],
                    "pair_count": r["Pair_Count"],
                    "cleaned_MCS_proxy": r["Mechanistic_Convergence_Score"],
                    "clinical_signal_proxy": r["Translational_Maturity_Index"],
                    "frontier_score_proxy": r["Gap_Score"],
                    "risk_reason": "|".join(reasons),
                }
            )
    warnings = pd.DataFrame(rows)
    warnings.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/low_count_warning_pairs_selected_death_modes_2000_2025.csv", index=False)
    reasons = warnings["risk_reason"].str.split("|").explode().value_counts().head(8) if not warnings.empty else pd.Series(dtype=int)
    modes = warnings["death_mode"].value_counts().reindex(DEATH_ORDER).dropna() if not warnings.empty else pd.Series(dtype=int)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    panel_label(axes[0], "A")
    axes[0].barh(reasons.index[::-1], reasons.values[::-1], color="#D28A2E")
    axes[0].set_xlabel("Flagged pairs")
    axes[0].set_title("Low-count inflation risk categories", loc="left", weight="bold")
    panel_label(axes[1], "B")
    axes[1].barh(modes.index[::-1], modes.values[::-1], color=[MODE_COLORS.get(m, "#888888") for m in modes.index[::-1]])
    axes[1].set_xlabel("Flagged pairs")
    axes[1].set_title("Warnings by concept", loc="left", weight="bold")
    for ax in axes:
        clean_axes(ax)
    savefig(fig, "Supplementary_Figure_S8_low_count_instability_selected_death_modes_2000_2025", False)


def s9(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"]
    rows = []
    for mode, sub in pair.groupby("Death_Mode"):
        total_pairs = len(sub)
        retained10 = (sub["Pair_Count"] >= 10).sum()
        retained_pmid_frac = sub[sub["Pair_Count"] >= 10]["Pair_Count"].sum() / sub["Pair_Count"].sum() if sub["Pair_Count"].sum() else 0
        low_frac = (sub["Pair_Count"] < 5).mean()
        mcs_drift = abs(sub.loc[sub["Pair_Count"] >= 10, "Mechanistic_Convergence_Score"].mean() - sub["Mechanistic_Convergence_Score"].mean())
        frontier_drift = abs(sub.loc[sub["Pair_Count"] >= 10, "Gap_Score"].mean() - sub["Gap_Score"].mean())
        retained_pair_fraction = retained10 / total_pairs if total_pairs else 0
        instability = (1 - retained_pair_fraction) + (1 - retained_pmid_frac) + low_frac + np.nan_to_num(mcs_drift) + np.nan_to_num(frontier_drift)
        rows.append(
            {
                "death_mode": mode,
                "retained_pair_fraction": retained_pair_fraction,
                "retained_pmid_fraction": retained_pmid_frac,
                "low_count_pair_fraction": low_frac,
                "MCS_drift": mcs_drift,
                "frontier_drift": frontier_drift,
                "instability_score": instability,
            }
        )
    df = pd.DataFrame(rows).sort_values("instability_score")
    df.to_csv(V5 / "04_supplementary_figures_regenerated/source_data/death_mode_sensitivity_summary_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), constrained_layout=True, gridspec_kw={"width_ratios": [0.9, 1.1]})
    panel_label(axes[0], "A")
    y = np.arange(len(df))
    axes[0].hlines(y, 0, df["instability_score"], color="#D0D0D0")
    axes[0].scatter(df["instability_score"], y, color=[MODE_COLORS.get(m, "#777777") for m in df["death_mode"]], s=42)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(df["death_mode"])
    axes[0].set_xlabel("Instability score")
    axes[0].set_title("Death-mode instability ranking", loc="left", weight="bold")
    clean_axes(axes[0])
    panel_label(axes[1], "B")
    cols = ["retained_pair_fraction", "retained_pmid_fraction", "low_count_pair_fraction", "MCS_drift", "frontier_drift"]
    sns.heatmap(df.set_index("death_mode")[cols], ax=axes[1], cmap="vlag", center=0, linewidths=0.25, cbar_kws={"shrink": 0.55})
    axes[1].set_title("Instability components", loc="left", weight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    savefig(fig, "Supplementary_Figure_S9_death_mode_instability_selected_death_modes_2000_2025", False)


def regenerate_figures(data: dict[str, pd.DataFrame]) -> None:
    figure1(data)
    figure2(data)
    figure3(data)
    figure4(data)
    figure5(data)
    s1(data)
    s2(data)
    s3(data)
    s4_generated = s4(data)
    s5(data)
    s6(data)
    s7(data)
    s8(data)
    s9(data)
    if not s4_generated:
        write_text(V5 / "04_supplementary_figures_regenerated/Supplementary_Figure_S4_missing_source_report.md", "# Supplementary Figure S4 Missing Source\n\nS4 was not regenerated because denominator-normalized source data were incomplete or unavailable after filtering.\n")


def build_manuscript_patch(data: dict[str, pd.DataFrame]) -> None:
    pair = data["pair"]
    tmi = data["tmi"]
    yearly = data["yearly"]
    stats = data["rank_stats"].iloc[0]
    s8_path = V5 / "04_supplementary_figures_regenerated/source_data/low_count_warning_pairs_selected_death_modes_2000_2025.csv"
    low_n = len(read_csv(s8_path)) if s8_path.exists() else 0
    s6_path = V5 / "04_supplementary_figures_regenerated/source_data/tmi_summary_selected_death_modes_2000_2025.csv"
    tmi_summary = read_csv(s6_path)
    text = f"""
    # Eight-Concept 2000-2025 Manuscript Patch

    This patch updates the figure-linked Results to the selected death-mode, {START_YEAR}-{END_YEAR} atlas.

    ## Scale sentence

    After applying the final eight-mode scope and excluding records outside {START_YEAR}-{END_YEAR}, the regenerated atlas contains {len(DEATH_ORDER)} RCD-associated concepts, {pair['Disease'].nunique():,} disease terms with at least one positive pair, {len(pair):,} positive death-mode disease pairs, {len(tmi):,} article-stage pair-PMID records, and {tmi['pmid'].nunique():,} unique PubMed records.

    ## Temporal sentence

    Figure 2 should now be described as a {START_YEAR}-{END_YEAR} terminology-adoption timeline. The retained chronology emphasizes inflammatory death, ferroptosis/redox expansion, and newer metal-stress or integrative terminology.

    ## Oncology-sensitivity sentence

    In the eight-concept atlas, oncology-sensitive filtering retained {int((~data['oncology']['oncology_flag']).sum()):,} non-oncology disease terms and flagged {int(data['oncology']['oncology_flag'].sum()):,} oncology/neoplasm terms. Global hub rankings remained correlated after neoplasm removal (Spearman = {stats['spearman']:.3f}), while top-hub identity overlap was limited (top25 Jaccard = {stats['top25_jaccard']:.3f}).

    ## Robustness sentence

    Low-count high-risk pairs after filtering: {low_n:,}. These pairs should remain exploratory and should not be presented as validated disease mechanisms.

    ## Figure 5 sentence

    Figure 5 remains a cleaned mechanism-vocabulary figure and should retain the sentence: These signals reflect literature-level mechanism vocabulary, not experimental proof of pathway equivalence.
    """
    write_text(V5 / "05_manuscript_patches/selected_death_modes_results_patch.md", text)

    legend = f"""
    # Eight-Concept Figure Legend Patch

    ## Figure 1
    Eight-concept PubMed atlas construction and disease-literature scale. The workflow begins with eight ontology-qualified death-mode query definitions and the cleaned MeSH-derived disease-term set, followed by pairwise PubMed retrieval, filtering to {START_YEAR}-{END_YEAR}, pair-PMID links, and audited atlas layers. Co-mention indicates literature association, not causal evidence.

    ## Figure 2
    Historical waves of RCD terminology adoption, {START_YEAR}-{END_YEAR}. Curves show annual unique PMID publication counts after 2000-2025 filtering. These are terminology-adoption waves, not biological emergence.

    ## Figure 3
    Disease-literature selectivity across eight RCD-associated concepts. Bubble area represents {START_YEAR}-{END_YEAR} pair count; color represents disease-literature selectivity.

    ## Figure 4
    Oncology-sensitive literature hubs in the eight-concept atlas. Hubs are literature-structure descriptors. Spearman indicates global rank similarity; top25 Jaccard indicates overlap of top hub identities.

    ## Figure 5
    Cleaned mechanism-vocabulary axes across eight-concept RCD-associated disease literature. Vocabulary-level signals; not experimental proof of pathway equivalence.
    """
    write_text(V5 / "05_manuscript_patches/selected_death_modes_figure_legend_patch.md", legend)


def file_inventory() -> None:
    cleanup_appledouble()
    rows = []
    for folder, placement in [
        ("03_main_figures_regenerated", "main"),
        ("04_supplementary_figures_regenerated", "supplementary"),
    ]:
        base = V5 / folder
        for pdf in sorted((base / "pdf").glob("*.pdf")):
            stem = pdf.stem
            rows.append(
                {
                    "figure_id": stem,
                    "placement": placement,
                    "pdf": str(pdf),
                    "png": str(base / "png" / f"{stem}.png"),
                    "tiff": str(base / "tiff" / f"{stem}.tiff"),
                    "pdf_exists": pdf.exists(),
                    "png_exists": (base / "png" / f"{stem}.png").exists(),
                    "tiff_exists": (base / "tiff" / f"{stem}.tiff").exists(),
                    "source_filter": f"selected death modes; {START_YEAR}-{END_YEAR}",
                }
            )
    write_csv(
        V5 / "10_final_qc/figure_inventory_selected_death_modes_2000_2025.csv",
        rows,
        ["figure_id", "placement", "pdf", "png", "tiff", "pdf_exists", "png_exists", "tiff_exists", "source_filter"],
    )


def markdown_text_check() -> list[str]:
    bad_year = []
    # Use text-level source and logs where available. Binary PDF text extraction is optional in QC.
    for p in list((V5 / "05_manuscript_patches").glob("*.md")) + list((V5 / "00_filter_contract").glob("*.md")):
        txt = p.read_text(encoding="utf-8", errors="ignore")
        if re.search(r"\b2026\b|199[0-9]\b", txt):
            bad_year.append(str(p))
    return bad_year


def build_qc(data: dict[str, pd.DataFrame]) -> None:
    cleanup_appledouble()
    file_inventory()
    bad_year = markdown_text_check()
    main_pdf = [p for p in (V5 / "03_main_figures_regenerated/pdf").glob("*.pdf") if not p.name.startswith("._")]
    supp_pdf = [p for p in (V5 / "04_supplementary_figures_regenerated/pdf").glob("*.pdf") if not p.name.startswith("._")]
    pair = data["pair"]
    tmi = data["tmi"]
    s4_exists = bool(list((V5 / "04_supplementary_figures_regenerated/pdf").glob("Supplementary_Figure_S4*.pdf")))
    text = f"""
    # Eight-Concept 2000-2025 Package QC Report

    Created: {TODAY}

    ## Filter

    - Scope: final eight selected regulated cell death modes
    - Year window: {START_YEAR}-{END_YEAR}
    - Death-mode concepts retained: {len(DEATH_ORDER)}
    - Disease terms retained: {pair['Disease'].nunique():,}
    - Positive pairs retained: {len(pair):,}
    - Pair-PMID stage rows retained: {len(tmi):,}
    - Unique PMIDs retained: {tmi['pmid'].nunique():,}

    ## Figures

    - Main figure PDFs: {len(main_pdf)}
    - Supplementary figure PDFs: {len(supp_pdf)}
    - Supplementary Figure S4 generated: {'yes' if s4_exists else 'no'}
    - All figures exported as PDF, PNG, and TIFF: {'yes' if all((p.with_suffix('.png').exists() for p in (V5/'03_main_figures_regenerated/pdf').glob('*.pdf'))) else 'see inventory'}

    ## Interpretation safety

    - Co-mention remains literature association only.
    - Cleaned mechanism vocabulary remains vocabulary-level only.
    - Evidence-stage audit is not translational maturity.
    - Hub scores are literature-structure descriptors.

    ## Text QC

    - Out-of-window year mentions in generated Markdown: {len(bad_year)}

    ## Notes

    Figure-level binary text should still be visually inspected before submission. If you decide to make the eight-concept filtered version the manuscript baseline, update abstract/results/methods counts using `05_manuscript_patches/selected_death_modes_results_patch.md`.
    """
    write_text(V5 / "10_final_qc/FINAL_SELECTED_DEATH_MODES_QC_REPORT.md", text)
    write_text(V5 / "10_final_qc/FINAL_FILE_TREE.txt", "\n".join(str(p.relative_to(V5)) for p in sorted(V5.rglob("*")) if p.is_file()))


def copy_script() -> None:
    dest = V5 / "02_scripts/build_v5_selected_death_modes_2000_2025_package.py"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__), dest)


def cleanup_appledouble() -> None:
    if not V5.exists():
        return
    for p in V5.rglob("._*"):
        if p.is_file():
            p.unlink()


def main() -> int:
    if not V4.exists():
        raise SystemExit(f"Missing v4 package: {V4}")
    mkdirs()
    data = compute_sources()
    save_sources(data)
    build_filter_contract(data)
    regenerate_figures(data)
    build_manuscript_patch(data)
    copy_script()
    build_qc(data)
    cleanup_appledouble()
    print(f"Built selected death-mode 2000-2025 package: {V5}")
    print(f"Main figures: {len([p for p in (V5 / '03_main_figures_regenerated/pdf').glob('*.pdf') if not p.name.startswith('._')])}")
    print(f"Supplementary figures: {len([p for p in (V5 / '04_supplementary_figures_regenerated/pdf').glob('*.pdf') if not p.name.startswith('._')])}")
    print(f"Positive pairs: {len(data['pair'])}")
    print(f"Unique PMIDs: {data['tmi']['pmid'].nunique()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
