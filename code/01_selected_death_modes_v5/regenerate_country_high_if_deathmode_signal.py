#!/usr/bin/env python3
"""Draw country-level high-JIF article signal split by death mode.

This figure uses the v5 selected death-mode article-stage table, enriched PubMed
affiliation countries, and CopyofImpactFactor2024.csv. It does not overwrite
the existing country/journal context figures.
"""

from __future__ import annotations

import os
import re
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/celldeath_mplconfig_v5_high_if_country")

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", PROJECT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
MAIN = V5 / "03_main_figures_regenerated"
ARTICLE_STAGE = V5 / "01_source_data_filtered" / "article_stage_records_selected_death_modes_2000_2025.csv"
ENRICHED_META = PROJECT / "results" / "article_metadata_emerging_article_nonapoptosis.csv"
JIF_TABLE = PROJECT / "CopyofImpactFactor2024.csv"

STEM = "Figure_country_high_if_death_mode_signal_selected_death_modes_2000_2025"
START_YEAR = 2000
END_YEAR = 2025
JIF_THRESHOLD = 10.0
PLOT_N = 20

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
        "axes.titlesize": 8.5,
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


def normalize_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).upper()
    text = text.replace("&", " AND ")
    text = re.sub(r"[\u2010-\u2015-]", " ", text)
    text = re.sub(r"[^A-Z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    text = re.sub(r"^THE\s+", "", text)
    return text


def journal_variants(value: object) -> list[str]:
    if pd.isna(value):
        return []
    raw = str(value).strip()
    if not raw:
        return []
    variants = [raw]
    variants.append(re.split(r"\s*=\s*", raw, maxsplit=1)[0])
    variants.append(re.split(r"\s*:\s*", raw, maxsplit=1)[0])
    variants.append(re.sub(r"\s*\([^)]*\)", "", raw))
    variants.append(re.split(r"\s*:\s*", re.sub(r"\s*\([^)]*\)", "", raw), maxsplit=1)[0])
    seen = []
    for item in variants:
        norm = normalize_text(item)
        if norm and norm not in seen:
            seen.append(norm)
    return seen


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


def split_list(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [p.strip() for p in re.split(r"\s*[|;]\s*", text) if p.strip()]


def compact_label(value: str, width: int = 34) -> str:
    value = str(value)
    if len(value) <= width:
        return value
    return textwrap.shorten(value, width=width, placeholder="...")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.075, 1.04, label, transform=ax.transAxes, weight="bold", fontsize=9, va="top")


def clean_axes(ax: plt.Axes, proportion: bool = False) -> None:
    ax.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax.set_axisbelow(True)
    if proportion:
        ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0, decimals=0))
    else:
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))


def savefig(fig: plt.Figure, stem: str) -> None:
    for fmt in ["pdf", "png", "tiff"]:
        out = MAIN / fmt / f"{stem}.{fmt}"
        out.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "tiff":
            fig.savefig(out, bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        else:
            fig.savefig(out, bbox_inches="tight", dpi=600)
    plt.close(fig)


def load_article_data() -> pd.DataFrame:
    article = pd.read_csv(ARTICLE_STAGE, dtype={"pmid": str}, low_memory=False)
    article["publication_year_num"] = pd.to_numeric(
        article.get("publication_year_num", article.get("publication_year")), errors="coerce"
    )
    article = article[
        article["death_mode"].isin(DEATH_ORDER)
        & article["publication_year_num"].between(START_YEAR, END_YEAR, inclusive="both")
    ].copy()
    article["pmid"] = article["pmid"].astype(str)

    meta = pd.read_csv(ENRICHED_META, dtype={"PMID": str}, low_memory=False)
    keep = [c for c in ["PMID", "Countries", "Journal"] if c in meta.columns]
    meta = meta[keep].drop_duplicates("PMID").rename(columns={"PMID": "pmid", "Journal": "meta_journal"})

    article = article.merge(meta, on="pmid", how="left")
    article["journal_for_if"] = article["journal"].where(article["journal"].notna(), article.get("meta_journal"))
    return article


def load_jif_lookup() -> pd.DataFrame:
    jif = pd.read_csv(JIF_TABLE, low_memory=False)
    required = {"Name", "Abbr Name", "JIF"}
    missing = sorted(required.difference(jif.columns))
    if missing:
        raise ValueError(f"{JIF_TABLE} missing required columns: {missing}")
    jif["JIF"] = pd.to_numeric(jif["JIF"], errors="coerce")
    jif = jif[jif["JIF"].notna()].copy()

    rows = []
    for _, row in jif.iterrows():
        for source_col, method in [("Name", "journal_name"), ("Abbr Name", "journal_abbreviation")]:
            norm = normalize_text(row.get(source_col, ""))
            if norm:
                rows.append(
                    {
                        "journal_key": norm,
                        "IF_Name": row.get("Name"),
                        "IF_Abbr": row.get("Abbr Name"),
                        "JIF": row.get("JIF"),
                        "JIF5Years": row.get("JIF5Years"),
                        "Category": row.get("Category"),
                        "match_method": method,
                    }
                )
    lookup = pd.DataFrame(rows)
    lookup = lookup.sort_values("JIF", ascending=False).drop_duplicates("journal_key", keep="first")
    return lookup


def match_journal_to_jif(article: pd.DataFrame, lookup: pd.DataFrame) -> pd.DataFrame:
    unique_journals = (
        article[["journal_for_if"]]
        .dropna()
        .drop_duplicates()
        .rename(columns={"journal_for_if": "Journal"})
        .reset_index(drop=True)
    )
    matches = []
    lookup_index = lookup.set_index("journal_key")
    for journal in unique_journals["Journal"]:
        variants = journal_variants(journal)
        hit = None
        hit_variant = ""
        for variant in variants:
            if variant in lookup_index.index:
                hit = lookup_index.loc[variant]
                hit_variant = variant
                break
        if hit is None:
            matches.append(
                {
                    "Journal": journal,
                    "Journal_norm": variants[0] if variants else "",
                    "matched_yes_no": "no",
                    "match_variant": "",
                    "IF_Name": "",
                    "IF_Abbr": "",
                    "JIF": np.nan,
                    "JIF5Years": np.nan,
                    "Category": "",
                    "match_method": "unmatched",
                }
            )
        else:
            matches.append(
                {
                    "Journal": journal,
                    "Journal_norm": variants[0] if variants else "",
                    "matched_yes_no": "yes",
                    "match_variant": hit_variant,
                    "IF_Name": hit["IF_Name"],
                    "IF_Abbr": hit["IF_Abbr"],
                    "JIF": hit["JIF"],
                    "JIF5Years": hit["JIF5Years"],
                    "Category": hit["Category"],
                    "match_method": hit["match_method"],
                }
            )
    matched = pd.DataFrame(matches)
    article = article.merge(matched, left_on="journal_for_if", right_on="Journal", how="left")
    article["high_if_yes_no"] = np.where(article["JIF"] >= JIF_THRESHOLD, "yes", "no")
    return article


def expand_country_records(article: pd.DataFrame) -> pd.DataFrame:
    rows = []
    high = article[article["high_if_yes_no"].eq("yes")].copy()
    base_cols = [
        "pmid",
        "death_mode",
        "disease_term",
        "publication_year_num",
        "journal_for_if",
        "IF_Name",
        "JIF",
        "JIF5Years",
        "Category",
        "match_method",
        "Countries",
    ]
    high = high[[c for c in base_cols if c in high.columns]].drop_duplicates()
    for row in high.itertuples(index=False):
        countries = sorted({normalize_country(c) for c in split_list(getattr(row, "Countries", ""))})
        countries = [c for c in countries if c and c.lower() != "missing"]
        for country in countries:
            rows.append(
                {
                    "pmid": row.pmid,
                    "death_mode": row.death_mode,
                    "disease_term": getattr(row, "disease_term", ""),
                    "publication_year": getattr(row, "publication_year_num", ""),
                    "journal": getattr(row, "journal_for_if", ""),
                    "IF_Name": getattr(row, "IF_Name", ""),
                    "JIF": getattr(row, "JIF", np.nan),
                    "JIF5Years": getattr(row, "JIF5Years", np.nan),
                    "Category": getattr(row, "Category", ""),
                    "match_method": getattr(row, "match_method", ""),
                    "Country": country,
                }
            )
    return pd.DataFrame(rows)


def build_country_tables(expanded: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if expanded.empty:
        return pd.DataFrame(), pd.DataFrame()
    counts = (
        expanded.drop_duplicates(["pmid", "Country", "death_mode"])
        .groupby(["Country", "death_mode"])["pmid"]
        .nunique()
        .unstack(fill_value=0)
    )
    for mode in DEATH_ORDER:
        if mode not in counts.columns:
            counts[mode] = 0
    counts = counts[DEATH_ORDER].astype(float)
    counts["high_if_article_mode_count"] = counts[DEATH_ORDER].sum(axis=1)
    unique_pmids = (
        expanded.drop_duplicates(["pmid", "Country"])
        .groupby("Country")["pmid"]
        .nunique()
        .rename("unique_high_if_pmids")
    )
    counts = counts.join(unique_pmids, how="left").fillna({"unique_high_if_pmids": 0})
    counts = counts.sort_values("high_if_article_mode_count", ascending=False)
    counts["rank"] = np.arange(1, len(counts) + 1)
    wide = counts.reset_index()

    long = wide.melt(
        id_vars=["Country", "high_if_article_mode_count", "unique_high_if_pmids", "rank"],
        value_vars=DEATH_ORDER,
        var_name="death_mode",
        value_name="high_if_pmid_count",
    )
    long["death_mode_proportion"] = np.where(
        long["high_if_article_mode_count"] > 0,
        long["high_if_pmid_count"] / long["high_if_article_mode_count"],
        0.0,
    )
    return wide, long


def plot_count_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    plot_df = df.sort_values("high_if_article_mode_count", ascending=False).head(PLOT_N)
    plot_df = plot_df.sort_values("high_if_article_mode_count", ascending=True).reset_index(drop=True)
    y = np.arange(len(plot_df))
    left = np.zeros(len(plot_df), dtype=float)
    for mode in DEATH_ORDER:
        vals = plot_df[mode].to_numpy(dtype=float)
        ax.barh(
            y,
            vals,
            left=left,
            color=MODE_COLORS[mode],
            edgecolor="white",
            linewidth=0.25,
            height=0.76,
        )
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([compact_label(v, 32) for v in plot_df["Country"]])
    xmax = float(plot_df["high_if_article_mode_count"].max())
    ax.set_xlim(0, xmax * 1.18)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))
    for yy, total, pmids in zip(y, plot_df["high_if_article_mode_count"], plot_df["unique_high_if_pmids"]):
        ax.text(total + xmax * 0.018, yy, f"n={int(total):,}", va="center", ha="left", fontsize=5.6, color="#263238")
    ax.set_xlabel("Article-mode PMID count in journals with JIF >= 10")
    ax.set_title("High-impact country affiliation signal", loc="left", weight="bold")
    panel_label(ax, "A")
    clean_axes(ax)


def plot_proportion_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    plot_df = df.sort_values("high_if_article_mode_count", ascending=False).head(PLOT_N)
    plot_df = plot_df.sort_values("high_if_article_mode_count", ascending=True).reset_index(drop=True)
    denom = plot_df["high_if_article_mode_count"].replace(0, np.nan).to_numpy(dtype=float)
    y = np.arange(len(plot_df))
    left = np.zeros(len(plot_df), dtype=float)
    for mode in DEATH_ORDER:
        vals = plot_df[mode].to_numpy(dtype=float) / denom
        vals = np.nan_to_num(vals, nan=0.0)
        ax.barh(
            y,
            vals,
            left=left,
            color=MODE_COLORS[mode],
            edgecolor="white",
            linewidth=0.25,
            height=0.76,
        )
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([compact_label(v, 32) for v in plot_df["Country"]])
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.25))
    ax.set_xlabel("Death-mode proportion within JIF >= 10 country signal")
    ax.set_title("Death-mode composition", loc="left", weight="bold")
    panel_label(ax, "B")
    clean_axes(ax, proportion=True)


def write_outputs(article_if: pd.DataFrame, expanded: pd.DataFrame, wide: pd.DataFrame, long: pd.DataFrame) -> None:
    source_dir = MAIN / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)

    long.to_csv(source_dir / f"{STEM}_long_source_data.csv", index=False)
    wide.to_csv(source_dir / f"{STEM}_country_wide_source_data.csv", index=False)
    expanded.to_csv(source_dir / f"{STEM}_expanded_pmid_country_source_data.csv", index=False)

    audit_cols = [
        "pmid",
        "death_mode",
        "disease_term",
        "publication_year_num",
        "journal_for_if",
        "Countries",
        "matched_yes_no",
        "match_method",
        "IF_Name",
        "JIF",
        "JIF5Years",
        "Category",
        "high_if_yes_no",
    ]
    audit = article_if[[c for c in audit_cols if c in article_if.columns]].drop_duplicates()
    audit.to_csv(source_dir / f"{STEM}_article_jif_match_audit.csv", index=False)

    unmatched = (
        article_if[article_if["matched_yes_no"].ne("yes")]
        .groupby("journal_for_if", dropna=False)["pmid"]
        .nunique()
        .reset_index(name="unique_pmids")
        .sort_values("unique_pmids", ascending=False)
    )
    unmatched.to_csv(source_dir / f"{STEM}_unmatched_journals.csv", index=False)

    total_records = len(article_if.drop_duplicates(["pmid", "death_mode", "journal_for_if"]))
    matched_records = len(article_if[article_if["matched_yes_no"].eq("yes")].drop_duplicates(["pmid", "death_mode", "journal_for_if"]))
    high_if_records = len(article_if[article_if["high_if_yes_no"].eq("yes")].drop_duplicates(["pmid", "death_mode", "journal_for_if"]))
    qc = pd.DataFrame(
        [
            {"metric": "article_stage_source", "value": str(ARTICLE_STAGE)},
            {"metric": "enriched_metadata_source", "value": str(ENRICHED_META)},
            {"metric": "jif_source", "value": str(JIF_TABLE)},
            {"metric": "jif_threshold", "value": f">= {JIF_THRESHOLD:g}"},
            {"metric": "unique_pmid_death_journal_records", "value": total_records},
            {"metric": "jif_matched_pmid_death_journal_records", "value": matched_records},
            {"metric": "jif_match_fraction", "value": matched_records / total_records if total_records else 0},
            {"metric": "high_if_pmid_death_journal_records", "value": high_if_records},
            {"metric": "country_expanded_high_if_records", "value": len(expanded)},
            {"metric": "countries_with_high_if_signal", "value": wide["Country"].nunique() if not wide.empty else 0},
            {"metric": "countries_plotted", "value": min(PLOT_N, len(wide))},
            {"metric": "count_unit", "value": "PMID-death-mode country full-counting signal"},
            {"metric": "country_policy", "value": "Affiliation-derived; Hong Kong, Macau, and Taiwan merged into China"},
            {"metric": "palette_source", "value": "Figure 2 MODE_COLORS"},
        ]
    )
    qc.to_csv(source_dir / f"{STEM}_QC.csv", index=False)


def main() -> int:
    article = load_article_data()
    lookup = load_jif_lookup()
    article_if = match_journal_to_jif(article, lookup)
    expanded = expand_country_records(article_if)
    wide, long = build_country_tables(expanded)
    if wide.empty:
        raise RuntimeError("No country-level records remained after JIF >= 10 filtering.")

    write_outputs(article_if, expanded, wide, long)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 6.0), constrained_layout=False, sharey=False)
    plot_count_panel(axes[0], wide)
    plot_proportion_panel(axes[1], wide)

    handles = [plt.Rectangle((0, 0), 1, 1, color=MODE_COLORS[m]) for m in DEATH_ORDER]
    fig.legend(
        handles,
        DEATH_ORDER,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.91),
        frameon=False,
        title="Cell death mode",
    )
    fig.suptitle(
        "Country-level affiliation signal in high-impact journals",
        x=0.02,
        y=0.975,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    fig.text(
        0.02,
        0.030,
        "Source: v5 selected death-mode 2000-2025 atlas, PubMed affiliation countries, and CopyofImpactFactor2024.csv. "
        "Countries use full counting for multi-country papers; counts are PMID-death-mode signals, not national performance metrics.",
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#555555",
    )
    fig.subplots_adjust(left=0.17, right=0.98, top=0.77, bottom=0.12, wspace=0.30)
    savefig(fig, STEM)

    print(f"Wrote {MAIN / 'pdf' / f'{STEM}.pdf'}")
    print(f"Wrote {MAIN / 'png' / f'{STEM}.png'}")
    print(f"Wrote source data under {MAIN / 'source_data'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
