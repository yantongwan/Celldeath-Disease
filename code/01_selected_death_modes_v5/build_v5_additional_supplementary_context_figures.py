#!/usr/bin/env python3
"""Build additional selected death-mode 2000-2025 supplementary context figures.

The figures generated here are visual analogues of older country, institution,
journal, frontier, and evidence-stage maps, but they are rebuilt against the
v5 selected death-mode 2000-2025 atlas whenever source data allow it.

No PubMed queries are made. No metrics are invented. Country and institution
figures depend on the existing enriched PMID metadata and are marked as
bibliometric context.
"""

from __future__ import annotations

import csv
import itertools
import math
import os
import re
import textwrap
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_additional")

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
SUPP = V5 / "04_supplementary_figures_regenerated"
SRC = SUPP / "source_data"
LOG_DIR = SUPP / "additional_context_logs"
SCRIPT_COPY = V5 / "02_scripts" / "build_v5_additional_supplementary_context_figures.py"

ARTICLE_STAGE = V5 / "01_source_data_filtered" / "article_stage_records_selected_death_modes_2000_2025.csv"
PAIR_METRICS = V5 / "01_source_data_filtered" / "pair_metrics_selected_death_modes_2000_2025.csv"
ENRICHED_META = PROJECT / "results" / "article_metadata_emerging_article_nonapoptosis.csv"

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

MUTED = {
    "blue": "#6C9EAE",
    "orange": "#D28A2E",
    "red": "#B9553C",
    "green": "#5C8F7B",
    "purple": "#6B6BAE",
    "gray": "#888888",
    "light_gray": "#D7D7D7",
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


def mkdirs() -> None:
    for d in ["pdf", "png", "tiff", "source_data", "additional_context_logs"]:
        (SUPP / d).mkdir(parents=True, exist_ok=True)
    SCRIPT_COPY.parent.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def savefig(fig: plt.Figure, stem: str) -> list[Path]:
    outputs = []
    for fmt in ["pdf", "png", "tiff"]:
        out = SUPP / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
        outputs.append(out)
    plt.close(fig)
    return outputs


def clean_axes(ax) -> None:
    ax.grid(axis="x", color="#E7E7E7", lw=0.45)
    ax.tick_params(length=2.5)


def panel_label(ax, label: str) -> None:
    ax.text(-0.06, 1.04, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def add_note(ax, note: str) -> None:
    ax.text(0.0, -0.18, note, transform=ax.transAxes, fontsize=5.8, color="#555555", va="top")


def zscore(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    sd = vals.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(vals)), index=series.index)
    return (vals - vals.mean()) / sd


def split_list(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = re.split(r"\s*[|;]\s*", text)
    return [p.strip() for p in parts if p.strip()]


def normalize_country(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value).strip())
    replacements = {
        "Hong Kong": "China",
        "Hong Kong SAR": "China",
        "Macao": "China",
        "Macau": "China",
        "Taiwan": "China",
        "Taiwan, Province of China": "China",
        "United States of America": "United States",
        "USA": "United States",
        "UK": "United Kingdom",
        "England": "United Kingdom",
        "Scotland": "United Kingdom",
        "Wales": "United Kingdom",
    }
    return replacements.get(text, text)


def clean_institution(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value).strip())
    text = re.sub(r"\s*\(.*?\)\s*", " ", text).strip()
    text = re.sub(r"\.$", "", text)
    generic_prefixes = (
        "department of ",
        "departments of ",
        "division of ",
        "laboratory of ",
        "key laboratory of ",
        "school of ",
        "faculty of ",
        "college of ",
    )
    if text.lower().startswith(generic_prefixes) and "," not in text:
        return ""
    if len(text) < 4:
        return ""
    return text


def compact_label(text: str, width: int = 34) -> str:
    text = str(text)
    if len(text) <= width:
        return text
    return text[: width - 1] + "…"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    article = pd.read_csv(ARTICLE_STAGE, dtype={"pmid": str}, low_memory=False)
    article["publication_year_num"] = pd.to_numeric(
        article.get("publication_year_num", article.get("publication_year")), errors="coerce"
    )
    article = article[
        article["death_mode"].isin(DEATH_ORDER)
        & article["publication_year_num"].between(START_YEAR, END_YEAR, inclusive="both")
    ].copy()
    article["pmid"] = article["pmid"].astype(str)

    pair = pd.read_csv(PAIR_METRICS, low_memory=False)
    pair = pair[pair["Death_Mode"].isin(DEATH_ORDER)].copy()

    meta = pd.read_csv(ENRICHED_META, dtype={"PMID": str}, low_memory=False)
    keep = [c for c in ["PMID", "PubYear", "Countries", "Institutions", "JournalCountry"] if c in meta.columns]
    meta = meta[keep].drop_duplicates("PMID")
    meta = meta.rename(columns={"PMID": "pmid"})
    article = article.merge(meta, on="pmid", how="left")
    return article, pair, meta


def expand_entity(article: pd.DataFrame, source_col: str, entity_col: str) -> pd.DataFrame:
    rows = []
    base = article[["pmid", "death_mode", source_col]].drop_duplicates()
    for row in base.itertuples(index=False):
        vals = split_list(getattr(row, source_col))
        if entity_col == "Country":
            vals = [normalize_country(v) for v in vals]
        elif entity_col == "Institution":
            vals = [clean_institution(v) for v in vals]
        vals = sorted({v for v in vals if v and v.lower() != "missing"})
        for val in vals:
            rows.append({"pmid": row.pmid, "death_mode": row.death_mode, entity_col: val})
    return pd.DataFrame(rows)


def journal_entity(article: pd.DataFrame) -> pd.DataFrame:
    df = article[["pmid", "death_mode", "journal"]].drop_duplicates().copy()
    df = df[df["journal"].notna()]
    df["Journal"] = df["journal"].astype(str).str.strip()
    df = df[df["Journal"].ne("")]
    return df[["pmid", "death_mode", "Journal"]]


def compute_rca(long_df: pd.DataFrame, entity_col: str, top_n: int, min_articles: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    if long_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    counts = (
        long_df.drop_duplicates(["pmid", "death_mode", entity_col])
        .groupby([entity_col, "death_mode"])["pmid"]
        .nunique()
        .unstack(fill_value=0)
    )
    for mode in DEATH_ORDER:
        if mode not in counts.columns:
            counts[mode] = 0
    counts = counts[DEATH_ORDER]
    totals = counts.sum(axis=1)
    counts = counts[totals >= min_articles]
    totals = counts.sum(axis=1)
    keep = totals.sort_values(ascending=False).head(top_n).index
    counts = counts.loc[keep]
    if counts.empty:
        return counts.reset_index(), counts.reset_index()
    entity_total = counts.sum(axis=1)
    mode_total = counts.sum(axis=0)
    grand = float(mode_total.sum())
    pseudo = 0.5
    expected = (mode_total + pseudo) / (grand + pseudo * len(DEATH_ORDER))
    rca = counts.copy().astype(float)
    for mode in DEATH_ORDER:
        entity_share = (counts[mode] + pseudo) / (entity_total + pseudo * len(DEATH_ORDER))
        rca[mode] = np.log2(entity_share / expected[mode])
    rca = rca.clip(lower=-4, upper=4)
    counts_out = counts.reset_index()
    rca_out = rca.reset_index()
    return counts_out, rca_out


def plot_rca_heatmap(rca: pd.DataFrame, entity_col: str, stem: str, title: str, note: str) -> bool:
    if rca.empty or entity_col not in rca.columns:
        return False
    mat = rca.set_index(entity_col)[DEATH_ORDER]
    mat.index = [compact_label(x, 46 if entity_col == "Institution" else 34) for x in mat.index]
    height = max(3.6, 0.22 * len(mat) + 1.2)
    width = 8.8 if entity_col == "Institution" else 7.4
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=False)
    sns.heatmap(
        mat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-3,
        vmax=3,
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "log2 relative citation/article association", "shrink": 0.55},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, loc="left", weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    add_note(ax, note)
    if entity_col == "Institution":
        fig.subplots_adjust(left=0.42, right=0.88, bottom=0.22, top=0.92)
    elif entity_col == "Journal":
        fig.subplots_adjust(left=0.34, right=0.88, bottom=0.22, top=0.92)
    else:
        fig.subplots_adjust(left=0.20, right=0.88, bottom=0.20, top=0.92)
    savefig(fig, stem)
    return True


def country_rca(article: pd.DataFrame, inventory: list[dict]) -> None:
    long = expand_entity(article, "Countries", "Country")
    counts, rca = compute_rca(long, "Country", top_n=22, min_articles=10)
    counts.to_csv(SRC / "additional_country_death_mode_counts_selected_death_modes_2000_2025.csv", index=False)
    rca.to_csv(SRC / "additional_country_death_mode_log2_rca_selected_death_modes_2000_2025.csv", index=False)
    ok = plot_rca_heatmap(
        rca,
        "Country",
        "Additional_Figure_country_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "Country-level RCD literature associations in the final eight-mode atlas",
        "Affiliation-derived country context; Hong Kong, Macau, and Taiwan are merged into China.",
    )
    add_inventory(inventory, "country_death_mode_log2_rca_heatmap", ok, "Computed from enriched PMID affiliations merged to v5 article-stage records.")


def institution_rca_and_productivity(article: pd.DataFrame, inventory: list[dict]) -> None:
    long = expand_entity(article, "Institutions", "Institution")
    counts, rca = compute_rca(long, "Institution", top_n=28, min_articles=15)
    counts.to_csv(SRC / "additional_institution_death_mode_counts_selected_death_modes_2000_2025.csv", index=False)
    rca.to_csv(SRC / "additional_institution_death_mode_log2_rca_selected_death_modes_2000_2025.csv", index=False)
    ok_heat = plot_rca_heatmap(
        rca,
        "Institution",
        "Additional_Figure_institution_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "Institution-level RCD literature associations in the final eight-mode atlas",
        "Institution context is affiliation-derived and should not be interpreted as institutional quality.",
    )
    add_inventory(inventory, "institution_death_mode_log2_rca_heatmap", ok_heat, "Computed from enriched PMID institution metadata merged to v5 article-stage records.")

    if long.empty:
        add_inventory(inventory, "institution_article_productivity", False, "No institution metadata available.")
        return
    summary = (
        long.drop_duplicates(["pmid", "Institution"])
        .groupby("Institution", as_index=False)
        .agg(article_n=("pmid", "nunique"))
        .merge(
            long.drop_duplicates(["pmid", "death_mode", "Institution"])
            .groupby("Institution")["death_mode"]
            .nunique()
            .rename("death_mode_breadth"),
            on="Institution",
            how="left",
        )
        .sort_values("article_n", ascending=False)
    )
    summary.to_csv(SRC / "additional_institution_article_productivity_selected_death_modes_2000_2025.csv", index=False)
    plot = summary.head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.6, 5.8), constrained_layout=True)
    colors = plt.cm.viridis((plot["death_mode_breadth"] - plot["death_mode_breadth"].min()) / max(1, plot["death_mode_breadth"].max() - plot["death_mode_breadth"].min()))
    ax.barh(plot["Institution"].map(lambda x: compact_label(x, 42)), plot["article_n"], color=colors)
    ax.set_xlabel("Unique atlas PMIDs")
    ax.set_title("Institutional publication context in the eight-concept atlas", loc="left", weight="bold")
    clean_axes(ax)
    add_note(ax, "Affiliation-derived context only; not an institution quality or performance metric.")
    savefig(fig, "Additional_Figure_institution_article_productivity_selected_death_modes_2000_2025")
    add_inventory(inventory, "institution_article_productivity", True, "Top affiliation-derived institutions by unique filtered atlas PMIDs.")


def institution_network(article: pd.DataFrame, inventory: list[dict]) -> None:
    long = expand_entity(article, "Institutions", "Institution")
    if long.empty:
        add_inventory(inventory, "institution_collaboration_network", False, "No institution metadata available.")
        return
    top_nodes = (
        long.drop_duplicates(["pmid", "Institution"]).groupby("Institution")["pmid"].nunique().sort_values(ascending=False).head(30)
    )
    top_set = set(top_nodes.index)
    pmid_to_inst = (
        long[long["Institution"].isin(top_set)]
        .drop_duplicates(["pmid", "Institution"])
        .groupby("pmid")["Institution"]
        .apply(lambda s: sorted(set(s)))
    )
    edge_counter: Counter[tuple[str, str]] = Counter()
    for insts in pmid_to_inst:
        if len(insts) < 2:
            continue
        for a, b in itertools.combinations(insts[:8], 2):
            edge_counter[tuple(sorted((a, b)))] += 1
    edges = pd.DataFrame(
        [{"Source": a, "Target": b, "Weight": w} for (a, b), w in edge_counter.items()]
    ).sort_values("Weight", ascending=False)
    edges = edges.head(55)
    edges.to_csv(SRC / "additional_institution_collaboration_edges_selected_death_modes_2000_2025.csv", index=False)
    top_nodes.rename("article_n").reset_index().rename(columns={"index": "Institution"}).to_csv(
        SRC / "additional_institution_collaboration_nodes_selected_death_modes_2000_2025.csv", index=False
    )
    if edges.empty:
        add_inventory(inventory, "institution_collaboration_network", False, "No multi-institution edges among top institutions.")
        return
    graph = nx.Graph()
    for inst, n in top_nodes.items():
        graph.add_node(inst, article_n=int(n))
    for row in edges.itertuples(index=False):
        graph.add_edge(row.Source, row.Target, weight=float(row.Weight))
    graph.remove_nodes_from([node for node, degree in dict(graph.degree()).items() if degree == 0])
    pos = nx.spring_layout(graph, seed=20260513, k=1.05, weight="weight", iterations=500, scale=1.0)
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    weights = np.array([d["weight"] for _, _, d in graph.edges(data=True)])
    if len(weights):
        widths = 0.3 + 2.2 * (weights / weights.max())
    else:
        widths = []
    nx.draw_networkx_edges(graph, pos, ax=ax, width=widths, edge_color="#BBBBBB", alpha=0.45)
    node_sizes = [80 + 11 * math.sqrt(graph.nodes[n]["article_n"]) for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=node_sizes, node_color="#6C9EAE", edgecolors="white", linewidths=0.6)
    label_nodes = set(top_nodes.loc[[n for n in top_nodes.index if n in graph.nodes]].head(12).index)
    for n in label_nodes:
        x, y = pos[n]
        ax.text(
            x,
            y + 0.035,
            compact_label(n, 24),
            fontsize=5.4,
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.72},
        )
    ax.set_title("Institutional co-affiliation network in the filtered atlas", loc="left", weight="bold")
    ax.margins(0.12)
    ax.set_axis_off()
    ax.text(0.0, -0.04, "Edges indicate co-occurrence in author affiliations; this is not a collaboration-quality metric.", transform=ax.transAxes, fontsize=5.8, color="#555555")
    savefig(fig, "Additional_Figure_institution_collaboration_network_selected_death_modes_2000_2025")
    add_inventory(inventory, "institution_collaboration_network", True, "Co-affiliation network among top institutions from enriched PMID metadata.")


def journal_context(article: pd.DataFrame, inventory: list[dict]) -> None:
    long = journal_entity(article)
    counts, rca = compute_rca(long, "Journal", top_n=28, min_articles=20)
    counts.to_csv(SRC / "additional_journal_death_mode_counts_selected_death_modes_2000_2025.csv", index=False)
    rca.to_csv(SRC / "additional_journal_death_mode_log2_rca_selected_death_modes_2000_2025.csv", index=False)
    ok = plot_rca_heatmap(
        rca,
        "Journal",
        "Additional_Figure_journal_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "Journal-level RCD literature associations in the final eight-mode atlas",
        "Atlas-internal journal context only; not a journal quality or impact metric.",
    )
    add_inventory(inventory, "journal_death_mode_log2_rca_heatmap", ok, "Computed directly from v5 article-stage journal metadata.")

    recent = article[article["publication_year_num"].between(2021, END_YEAR, inclusive="both")]
    j = recent[["pmid", "death_mode", "journal"]].drop_duplicates()
    j = j[j["journal"].notna()].copy()
    j["Journal"] = j["journal"].astype(str).str.strip()
    total = j.drop_duplicates(["pmid", "Journal"]).groupby("Journal")["pmid"].nunique().sort_values(ascending=False)
    keep = total.head(28).index
    mat_counts = (
        j[j["Journal"].isin(keep)]
        .groupby(["Journal", "death_mode"])["pmid"]
        .nunique()
        .unstack(fill_value=0)
    )
    for mode in DEATH_ORDER:
        if mode not in mat_counts.columns:
            mat_counts[mode] = 0
    mat_counts = mat_counts[DEATH_ORDER]
    share = mat_counts.div(mat_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    share = share.loc[total.loc[keep].index]
    share.reset_index().to_csv(SRC / "additional_journal_recent_death_mode_share_selected_death_modes_2021_2025.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    sns.heatmap(
        share,
        ax=ax,
        cmap="YlGnBu",
        vmin=0,
        vmax=max(0.25, float(np.nanpercentile(share.values, 95))),
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "Atlas-internal recent share", "shrink": 0.55},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Recent journal distribution of RCD-associated atlas articles", loc="left", weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    add_note(ax, "Recent window: 2021-2025. Shares are atlas-internal and denominator-dependent.")
    savefig(fig, "Additional_Figure_journal_recent_death_mode_share_heatmap_selected_death_modes_2021_2025")
    add_inventory(inventory, "journal_recent_death_mode_share_heatmap", True, "Recent atlas-internal death-mode share by journal from v5 article-stage records.")


def frontier_maps(pair: pd.DataFrame, inventory: list[dict]) -> None:
    df = pair.copy()
    for col in ["Pair_Count", "Gap_Score", "Death_Selectivity_Score", "Mechanistic_Convergence_Score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["frontier_proxy"] = zscore(df["Gap_Score"]) + zscore(df["Death_Selectivity_Score"]) - 0.20 * zscore(np.log10(df["Pair_Count"].fillna(0) + 1))
    df["label"] = df["Death_Mode"].astype(str) + " | " + df["Disease"].astype(str)
    df.to_csv(SRC / "additional_frontier_gap_scores_selected_death_modes_2000_2025.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for mode, sub in df.groupby("Death_Mode", sort=False):
        ax.scatter(
            sub["Death_Selectivity_Score"],
            sub["Gap_Score"],
            s=np.clip(np.sqrt(sub["Pair_Count"]) * 8, 10, 120),
            color=MODE_COLORS.get(mode, "#888888"),
            alpha=0.45,
            edgecolor="white",
            linewidth=0.25,
            label=mode,
        )
    top = df[df["Pair_Count"].ge(3)].nlargest(12, "frontier_proxy")
    for row in top.itertuples(index=False):
        ax.text(row.Death_Selectivity_Score, row.Gap_Score, compact_label(row.Disease, 24), fontsize=5.2, color="#333333")
    ax.axhline(0, color="#BBBBBB", lw=0.6)
    ax.axvline(0, color="#BBBBBB", lw=0.6)
    ax.set_xlabel("Disease-literature selectivity")
    ax.set_ylabel("Literature gap score")
    ax.set_title("Exploratory literature-gap landscape across all disease pairs", loc="left", weight="bold")
    ax.legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    clean_axes(ax)
    add_note(ax, "Point size scales with pair count. Low-count high-gap pairs are hypothesis-generating only.")
    savefig(fig, "Additional_Figure_frontier_gap_selected_death_modes_2000_2025")
    add_inventory(inventory, "frontier_gap", True, "Computed from v5 pair metrics using gap and disease-literature selectivity.")

    balanced = (
        df[df["Pair_Count"].ge(3)]
        .sort_values(["Death_Mode", "frontier_proxy"], ascending=[True, False])
        .groupby("Death_Mode", group_keys=False)
        .head(5)
        .copy()
    )
    balanced["display"] = balanced["Disease"].map(lambda x: compact_label(x, 34))
    balanced.to_csv(SRC / "additional_frontier_gap_balanced_top_pairs_selected_death_modes_2000_2025.csv", index=False)
    balanced = balanced.sort_values(["Death_Mode", "frontier_proxy"])
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    y = np.arange(len(balanced))
    colors = [MODE_COLORS.get(x, "#888888") for x in balanced["Death_Mode"]]
    ax.scatter(balanced["frontier_proxy"], y, s=np.clip(np.sqrt(balanced["Pair_Count"]) * 14, 20, 180), c=colors, edgecolor="white", linewidth=0.35)
    ax.set_yticks(y)
    ax.set_yticklabels(balanced["Death_Mode"].astype(str) + " | " + balanced["display"].astype(str), fontsize=5.6)
    ax.set_xlabel("Balanced exploratory literature-gap score")
    ax.set_title("Balanced top literature-gap pairs within each RCD concept", loc="left", weight="bold")
    clean_axes(ax)
    add_note(ax, "Restricted to Pair_Count >= 3 and top five pairs per death-mode concept.")
    savefig(fig, "Additional_Figure_frontier_gap_balanced_selected_death_modes_2000_2025")
    add_inventory(inventory, "frontier_gap_balanced", True, "Within-mode balanced top literature-gap pairs from v5 pair metrics.")


def evidence_stage_maps(article: pd.DataFrame, inventory: list[dict]) -> None:
    df = article.copy()
    df["audited_tmi_stage_score"] = pd.to_numeric(df.get("audited_tmi_stage_score"), errors="coerce")
    df["clinical_vocabulary_signal"] = df.get("clinical_vocabulary_signal_yes_no", "no").astype(str).str.lower().eq("yes")
    df["trial_supported"] = df.get("interventional_supported_by_publication_type_yes_no", "no").astype(str).str.lower().eq("yes")
    pair_stage = (
        df.groupby(["death_mode", "disease_term"], as_index=False)
        .agg(
            pair_count=("pmid", "nunique"),
            audited_stage_score_max=("audited_tmi_stage_score", "max"),
            audited_stage_score_mean=("audited_tmi_stage_score", "mean"),
            clinical_vocabulary_article_n=("clinical_vocabulary_signal", "sum"),
            interventional_supported_article_n=("trial_supported", "sum"),
        )
    )
    pair_stage["clinical_vocabulary_fraction"] = pair_stage["clinical_vocabulary_article_n"] / pair_stage["pair_count"].replace(0, np.nan)
    pair_stage["evidence_stage_signal"] = (
        pair_stage["audited_stage_score_max"].fillna(0)
        + 0.25 * np.log1p(pair_stage["interventional_supported_article_n"])
        + 0.15 * pair_stage["clinical_vocabulary_fraction"].fillna(0)
    )
    pair_stage["label"] = pair_stage["death_mode"] + " | " + pair_stage["disease_term"]
    pair_stage.to_csv(SRC / "additional_evidence_stage_pair_signal_selected_death_modes_2000_2025.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for mode, sub in pair_stage.groupby("death_mode", sort=False):
        ax.scatter(
            np.log10(sub["pair_count"] + 1),
            sub["clinical_vocabulary_fraction"],
            s=np.clip(np.sqrt(sub["interventional_supported_article_n"] + 1) * 12, 12, 120),
            color=MODE_COLORS.get(mode, "#888888"),
            alpha=0.48,
            edgecolor="white",
            linewidth=0.25,
            label=mode,
        )
    top = pair_stage.nlargest(12, "evidence_stage_signal")
    for row in top.itertuples(index=False):
        ax.text(math.log10(row.pair_count + 1), row.clinical_vocabulary_fraction, compact_label(row.disease_term, 22), fontsize=5.2)
    ax.set_xlabel("log10(pair PMIDs + 1)")
    ax.set_ylabel("Clinical-vocabulary article fraction")
    ax.set_title("Evidence-stage and clinical-vocabulary signal map", loc="left", weight="bold")
    ax.legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    clean_axes(ax)
    add_note(ax, "Audited evidence-stage signal only; not clinical-readiness evidence.")
    savefig(fig, "Additional_Figure_evidence_stage_signal_map_selected_death_modes_2000_2025")
    add_inventory(inventory, "evidence_stage_signal_map", True, "Audited evidence-stage replacement for the legacy clinical-facing map.")

    balanced = (
        pair_stage[pair_stage["pair_count"].ge(3)]
        .sort_values(["death_mode", "evidence_stage_signal"], ascending=[True, False])
        .groupby("death_mode", group_keys=False)
        .head(5)
        .copy()
    )
    balanced = balanced.sort_values(["death_mode", "evidence_stage_signal"])
    balanced.to_csv(SRC / "additional_evidence_stage_signal_balanced_top_pairs_selected_death_modes_2000_2025.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    y = np.arange(len(balanced))
    ax.scatter(
        balanced["evidence_stage_signal"],
        y,
        s=np.clip(np.sqrt(balanced["pair_count"]) * 11, 20, 160),
        c=[MODE_COLORS.get(x, "#888888") for x in balanced["death_mode"]],
        edgecolor="white",
        linewidth=0.35,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(balanced["death_mode"] + " | " + balanced["disease_term"].map(lambda x: compact_label(x, 34)), fontsize=5.6)
    ax.set_xlabel("Balanced evidence-stage signal")
    ax.set_title("Balanced evidence-stage signals within each RCD concept", loc="left", weight="bold")
    clean_axes(ax)
    add_note(ax, "Balanced view uses Pair_Count >= 3 and top five pairs per mode.")
    savefig(fig, "Additional_Figure_evidence_stage_signal_map_balanced_selected_death_modes_2000_2025")
    add_inventory(inventory, "evidence_stage_signal_map_balanced", True, "Balanced audited evidence-stage replacement for the legacy clinical-facing map.")


def rank_stability_summary(inventory: list[dict]) -> None:
    hub_path = SRC / "hub_rank_stability_selected_death_modes_2000_2025.csv"
    metric_path = SRC / "metric_rank_stability_selected_death_modes_2000_2025.csv"
    hub = pd.read_csv(hub_path) if hub_path.exists() else pd.DataFrame()
    metric = pd.read_csv(metric_path) if metric_path.exists() else pd.DataFrame()
    if hub.empty and metric.empty:
        add_inventory(inventory, "S2_hub_rank_stability_across_pair_count_thresholds", False, "Rank-stability source data missing.")
        return
    hub.to_csv(SRC / "additional_S2_hub_rank_stability_source_selected_death_modes_2000_2025.csv", index=False)
    metric.to_csv(SRC / "additional_S2_metric_rank_stability_source_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    if not hub.empty:
        panel_label(axes[0], "A")
        axes[0].plot(hub["comparison_threshold"], hub["spearman_correlation"], marker="o", color=MUTED["blue"])
        axes[0].plot(hub["comparison_threshold"], hub["top_25_jaccard"], marker="s", color=MUTED["orange"])
        axes[0].set_xticks(sorted(hub["comparison_threshold"].unique()))
        axes[0].set_ylim(0, 1.05)
        axes[0].set_xlabel("Pair_Count threshold")
        axes[0].set_ylabel("Stability versus >=1")
        axes[0].set_title("Hub-rank stability", loc="left", weight="bold")
        axes[0].legend(["Spearman", "Top-25 Jaccard"], frameon=False)
        clean_axes(axes[0])
    if not metric.empty:
        panel_label(axes[1], "B")
        piv = metric.pivot_table(index="metric_name", columns="comparison_threshold", values="spearman_correlation", aggfunc="mean")
        sns.heatmap(piv, ax=axes[1], cmap="YlGnBu", vmin=0, vmax=1, linewidths=0.25, cbar_kws={"label": "Spearman", "shrink": 0.55})
        axes[1].set_xlabel("Pair_Count threshold")
        axes[1].set_ylabel("")
        axes[1].set_title("Metric-rank stability", loc="left", weight="bold")
    fig.suptitle("Supplementary Figure S2 analogue. Rank stability across pair-count thresholds", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Additional_Figure_S2_hub_rank_stability_across_pair_count_thresholds_selected_death_modes_2000_2025")
    add_inventory(inventory, "S2_hub_rank_stability_across_pair_count_thresholds", True, "Regenerated from v5 robustness source data.")


def evidence_stage_summary(article: pd.DataFrame, inventory: list[dict]) -> None:
    original_guideline = int(article.get("original_stage", pd.Series(dtype=str)).astype(str).eq("guideline_evidence").sum())
    valid_guideline = int(article.get("valid_guideline_yes_no", pd.Series(dtype=str)).astype(str).str.lower().eq("yes").sum())
    guideline_ineligible = int(article.get("guideline_evidence_status", pd.Series(dtype=str)).astype(str).eq("guideline_ineligible_in_article_corpus").sum())
    trial_supported = int(article.get("interventional_supported_by_publication_type_yes_no", pd.Series(dtype=str)).astype(str).str.lower().eq("yes").sum())
    total_records = len(article)
    summary = pd.DataFrame(
        [
            {"category": "PMID-stage records audited", "article_stage_records": total_records},
            {"category": "Original guideline assignments", "article_stage_records": original_guideline},
            {"category": "Valid guideline evidence retained", "article_stage_records": valid_guideline},
            {"category": "Guideline-ineligible in article corpus", "article_stage_records": guideline_ineligible},
            {"category": "PublicationType-supported trials", "article_stage_records": trial_supported},
        ]
    )

    stage_order = [
        "basic_mechanism",
        "preclinical_animal",
        "human_observational",
        "interventional_trial",
        "guideline_ineligible_in_article_corpus",
    ]
    stage_labels = {
        "basic_mechanism": "Basic\nmechanism",
        "preclinical_animal": "Preclinical\nanimal",
        "human_observational": "Human\nobservational",
        "interventional_trial": "PublicationType-\nsupported trial",
        "guideline_ineligible_in_article_corpus": "Guideline-like\nineligible",
    }
    stage_colors = {
        "basic_mechanism": "#DDD1B6",
        "preclinical_animal": "#BFA47D",
        "human_observational": "#8FAFB8",
        "interventional_trial": "#3E6F98",
        "guideline_ineligible_in_article_corpus": "#DADDE2",
    }
    composition_count = (
        article.groupby(["death_mode", "audited_stage"], as_index=False)
        .agg(article_stage_records=("pmid", "size"))
    )
    total_by_mode = article.groupby("death_mode")["pmid"].size().rename("mode_stage_records")
    composition = composition_count.merge(total_by_mode, on="death_mode", how="left")
    composition["stage_share"] = composition["article_stage_records"] / composition["mode_stage_records"]
    complete_rows = []
    for mode in DEATH_ORDER:
        for stage in stage_order:
            sub = composition[(composition["death_mode"].eq(mode)) & (composition["audited_stage"].eq(stage))]
            if sub.empty:
                complete_rows.append(
                    {
                        "death_mode": mode,
                        "audited_stage": stage,
                        "article_stage_records": 0,
                        "mode_stage_records": int(total_by_mode.get(mode, 0)),
                        "stage_share": 0.0,
                    }
                )
            else:
                complete_rows.append(sub.iloc[0].to_dict())
    composition = pd.DataFrame(complete_rows)

    high_stage = (
        article.groupby("death_mode", as_index=False)
        .agg(
            mode_stage_records=("pmid", "size"),
            original_high_stage_text_signal=("original_high_stage_yes_no", lambda x: x.astype(str).str.lower().eq("yes").sum()),
            original_guideline_text_signal=("original_stage", lambda x: x.astype(str).eq("guideline_evidence").sum()),
            original_interventional_text_signal=("original_stage", lambda x: x.astype(str).eq("interventional_trial").sum()),
            publication_type_supported_trial=("interventional_supported_by_publication_type_yes_no", lambda x: x.astype(str).str.lower().eq("yes").sum()),
            valid_guideline_evidence=("valid_guideline_yes_no", lambda x: x.astype(str).str.lower().eq("yes").sum()),
        )
    )
    for col in [
        "original_high_stage_text_signal",
        "original_guideline_text_signal",
        "original_interventional_text_signal",
        "publication_type_supported_trial",
        "valid_guideline_evidence",
    ]:
        high_stage[col + "_share"] = high_stage[col] / high_stage["mode_stage_records"].replace(0, np.nan)
    high_stage["death_mode"] = pd.Categorical(high_stage["death_mode"], categories=DEATH_ORDER, ordered=True)
    high_stage = high_stage.sort_values("death_mode")

    summary.to_csv(SRC / "additional_S6_evidence_stage_audit_summary_selected_death_modes_2000_2025.csv", index=False)
    composition.to_csv(SRC / "additional_S6_evidence_stage_composition_selected_death_modes_2000_2025.csv", index=False)
    high_stage.to_csv(SRC / "additional_S6_high_stage_signal_after_audit_selected_death_modes_2000_2025.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.7), constrained_layout=False, gridspec_kw={"width_ratios": [1.45, 1.0]})
    short_mode_labels = {
        "Ferroptosis": "Ferroptosis",
        "Pyroptosis": "Pyroptosis",
        "NETosis": "NETosis",
        "Necroptosis": "Necroptosis",
        "Immunogenic cell death": "ICD",
        "Cuproptosis": "Cuproptosis",
        "PANoptosis": "PANoptosis",
        "Disulfidptosis": "Disulfidptosis",
    }
    y = np.arange(len(DEATH_ORDER))
    left = np.zeros(len(DEATH_ORDER))
    for stage in stage_order:
        vals = (
            composition[composition["audited_stage"].eq(stage)]
            .set_index("death_mode")
            .reindex(DEATH_ORDER)["stage_share"]
            .fillna(0)
            .to_numpy()
        )
        kwargs = {}
        if stage == "guideline_ineligible_in_article_corpus":
            kwargs = {"hatch": "///", "edgecolor": "white", "linewidth": 0.35}
        axes[0].barh(y, vals, left=left, height=0.78, color=stage_colors[stage], label=stage_labels[stage], **kwargs)
        left += vals
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([short_mode_labels.get(x, x) for x in DEATH_ORDER])
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Share of PMID-stage records")
    axes[0].set_title("A. Audited evidence-stage composition", loc="left", weight="bold")
    axes[0].legend(frameon=False, fontsize=5.6, loc="lower right", bbox_to_anchor=(1.0, 0.02))
    clean_axes(axes[0])

    hs = high_stage.set_index("death_mode").reindex(DEATH_ORDER)
    axes[1].barh(
        y,
        hs["original_high_stage_text_signal_share"],
        height=0.34,
        color="#D7DCE1",
        label="Original high-stage text signal",
    )
    axes[1].barh(
        y - 0.22,
        hs["publication_type_supported_trial_share"],
        height=0.18,
        color="#3E6F98",
        label="PublicationType-supported trial",
    )
    axes[1].set_yticks(y)
    axes[1].set_yticklabels([short_mode_labels.get(x, x) for x in DEATH_ORDER])
    axes[1].set_xlim(0, max(0.82, float(hs["original_high_stage_text_signal_share"].max()) * 1.05))
    axes[1].set_xlabel("Share")
    axes[1].set_title("B. High-stage signal after audit", loc="left", weight="bold")
    axes[1].legend(frameon=False, fontsize=5.7, loc="lower right")
    stats_text = (
        f"{total_records:,} PMID-stage records audited\n"
        f"{original_guideline:,} guideline assignments invalidated/ineligible\n"
        f"{valid_guideline:,} valid guideline evidence retained\n"
        f"{trial_supported:,} trial assignments supported by PublicationTypes"
    )
    axes[1].text(
        0.98,
        0.98,
        stats_text,
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color="#4A5668",
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#CBD5E1", "alpha": 0.92},
    )
    clean_axes(axes[1])
    fig.suptitle("Supplementary Figure S6 analogue. Evidence-stage audit and clinical-vocabulary signal", x=0.01, ha="left", fontsize=8.8, weight="bold")
    fig.text(0.01, 0.03, "Supplementary-only calibration: clinical-facing vocabulary is a literature signal; interventional evidence requires PubMed PublicationType support.", fontsize=5.9, color="#6B7C93")
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.20, top=0.82, wspace=0.32)
    savefig(fig, "Additional_Figure_S6_evidence_stage_audit_and_clinical_vocabulary_signal_selected_death_modes_2000_2025")
    add_inventory(inventory, "S6_evidence_stage_audit_and_clinical_vocabulary_signal", True, "Regenerated from v5 article-stage audit data.")


def death_mode_instability_summary(inventory: list[dict]) -> None:
    path = SRC / "death_mode_sensitivity_summary_selected_death_modes_2000_2025.csv"
    df = pd.read_csv(path) if path.exists() else pd.DataFrame()
    if df.empty:
        add_inventory(inventory, "S9_death_mode_instability_summary", False, "Death-mode sensitivity source data missing.")
        return
    df = df.copy()
    df["death_mode"] = pd.Categorical(df["death_mode"], categories=DEATH_ORDER, ordered=True)
    df = df.sort_values("instability_score", ascending=True)
    df.to_csv(SRC / "additional_S9_death_mode_instability_summary_selected_death_modes_2000_2025.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    panel_label(axes[0], "A")
    axes[0].barh(df["death_mode"], df["instability_score"], color=[MODE_COLORS.get(x, "#888888") for x in df["death_mode"].astype(str)])
    axes[0].set_xlabel("Instability score")
    axes[0].set_title("Death-mode instability ranking", loc="left", weight="bold")
    clean_axes(axes[0])
    panel_label(axes[1], "B")
    cols = ["retained_pair_fraction", "retained_pmid_fraction", "low_count_pair_fraction", "MCS_drift", "frontier_drift"]
    sns.heatmap(df.set_index("death_mode")[cols], ax=axes[1], cmap="vlag", center=0, linewidths=0.25, cbar_kws={"shrink": 0.55})
    axes[1].set_title("Instability components", loc="left", weight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    fig.suptitle("Supplementary Figure S9 analogue. Death-mode instability summary", x=0.01, ha="left", fontsize=10, weight="bold")
    savefig(fig, "Additional_Figure_S9_death_mode_instability_summary_selected_death_modes_2000_2025")
    add_inventory(inventory, "S9_death_mode_instability_summary", True, "Regenerated from v5 death-mode sensitivity data.")


def add_inventory(inventory: list[dict], figure_key: str, generated: bool, note: str) -> None:
    stem_map = {
        "country_death_mode_log2_rca_heatmap": "Additional_Figure_country_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "frontier_gap": "Additional_Figure_frontier_gap_selected_death_modes_2000_2025",
        "frontier_gap_balanced": "Additional_Figure_frontier_gap_balanced_selected_death_modes_2000_2025",
        "institution_article_productivity": "Additional_Figure_institution_article_productivity_selected_death_modes_2000_2025",
        "institution_collaboration_network": "Additional_Figure_institution_collaboration_network_selected_death_modes_2000_2025",
        "institution_death_mode_log2_rca_heatmap": "Additional_Figure_institution_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "journal_death_mode_log2_rca_heatmap": "Additional_Figure_journal_death_mode_log2_rca_heatmap_selected_death_modes_2000_2025",
        "journal_recent_death_mode_share_heatmap": "Additional_Figure_journal_recent_death_mode_share_heatmap_selected_death_modes_2021_2025",
        "evidence_stage_signal_map": "Additional_Figure_evidence_stage_signal_map_selected_death_modes_2000_2025",
        "evidence_stage_signal_map_balanced": "Additional_Figure_evidence_stage_signal_map_balanced_selected_death_modes_2000_2025",
        "S2_hub_rank_stability_across_pair_count_thresholds": "Additional_Figure_S2_hub_rank_stability_across_pair_count_thresholds_selected_death_modes_2000_2025",
        "S6_evidence_stage_audit_and_clinical_vocabulary_signal": "Additional_Figure_S6_evidence_stage_audit_and_clinical_vocabulary_signal_selected_death_modes_2000_2025",
        "S9_death_mode_instability_summary": "Additional_Figure_S9_death_mode_instability_summary_selected_death_modes_2000_2025",
    }
    stem = stem_map.get(figure_key, figure_key)
    inventory.append(
        {
            "figure_key": figure_key,
            "generated_yes_no": "yes" if generated else "no",
            "pdf": str(SUPP / "pdf" / f"{stem}.pdf") if generated else "MISSING",
            "png": str(SUPP / "png" / f"{stem}.png") if generated else "MISSING",
            "tiff": str(SUPP / "tiff" / f"{stem}.tiff") if generated else "MISSING",
            "note": note,
        }
    )


def write_inventory(inventory: list[dict]) -> None:
    fields = ["figure_key", "generated_yes_no", "pdf", "png", "tiff", "note"]
    out = LOG_DIR / "additional_supplementary_context_figure_inventory.csv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(inventory)
    generated = sum(1 for row in inventory if row["generated_yes_no"] == "yes")
    missing = [row["figure_key"] for row in inventory if row["generated_yes_no"] != "yes"]
    report = f"""# Additional Supplementary Context Figure Report

Created: {TODAY}

## Scope

- Package: v5 selected death-mode atlas
- Year window: {START_YEAR}-{END_YEAR}
- Retained death-mode concepts: {len(DEATH_ORDER)}
- Included death-mode concepts: {', '.join(DEATH_ORDER)}
- Requested figure analogues: {len(inventory)}
- Generated figure analogues: {generated}
- Missing or not generated: {len(missing)}

## Interpretation guardrails

- Country, institution, and collaboration figures are affiliation-derived bibliometric context only.
- Journal figures are atlas-internal journal context unless denominator-normalized source data are explicitly used.
- Frontier/gap figures are hypothesis-generating literature-structure maps, not validated biology.
- Evidence-stage figures replace legacy clinical-facing language and should not be read as clinical readiness.

## Missing figure keys

{chr(10).join('- ' + x for x in missing) if missing else 'None'}
"""
    write_text(LOG_DIR / "additional_supplementary_context_figure_report.md", report)


def qc_pdf_text() -> None:
    forbidden = [
        "translational maturity",
        "biological convergence",
        "validated clinical readiness",
        "proven mechanism",
        "true biological convergence",
    ]
    hits = []
    try:
        import subprocess

        for pdf in (SUPP / "pdf").glob("Additional_Figure*.pdf"):
            res = subprocess.run(["pdftotext", str(pdf), "-"], check=False, capture_output=True, text=True)
            txt = res.stdout.lower()
            for term in forbidden:
                if term in txt:
                    hits.append({"file": str(pdf), "term": term})
    except Exception as exc:
        hits.append({"file": "PDF_TEXT_CHECK_SKIPPED", "term": str(exc)})
    pd.DataFrame(hits).to_csv(LOG_DIR / "additional_supplementary_context_forbidden_text_check.csv", index=False)


def main() -> int:
    mkdirs()
    article, pair, _ = load_data()
    inventory: list[dict] = []

    country_rca(article, inventory)
    institution_rca_and_productivity(article, inventory)
    institution_network(article, inventory)
    journal_context(article, inventory)
    frontier_maps(pair, inventory)
    evidence_stage_maps(article, inventory)
    rank_stability_summary(inventory)
    evidence_stage_summary(article, inventory)
    death_mode_instability_summary(inventory)

    write_inventory(inventory)
    qc_pdf_text()
    if Path(__file__).resolve() != SCRIPT_COPY:
        SCRIPT_COPY.write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")

    generated = sum(1 for row in inventory if row["generated_yes_no"] == "yes")
    print(f"additional_figures_requested={len(inventory)}")
    print(f"additional_figures_generated={generated}")
    print(f"output_dir={SUPP}")
    print(f"inventory={LOG_DIR / 'additional_supplementary_context_figure_inventory.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
