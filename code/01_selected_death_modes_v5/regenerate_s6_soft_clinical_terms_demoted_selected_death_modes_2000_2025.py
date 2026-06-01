#!/usr/bin/env python3
"""Regenerate S6 with soft clinical-vocabulary terms demoted when appropriate.

This variant keeps the original S6 intact and writes a suffixed figure. The
rule change is deliberately conservative: serum/plasma/biomarker/survival
analysis only trigger reassignment from human_observational to
preclinical_animal when there is no stronger human-observational evidence and
there is explicit preclinical_animal text support.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", ROOT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
ARTICLE_STAGE = V5 / "01_source_data_filtered" / "article_stage_records_selected_death_modes_2000_2025.csv"
OUT_ROOT = V5 / "04_supplementary_figures_regenerated"
SRC = OUT_ROOT / "source_data"
STEM = "Additional_Figure_S6_evidence_stage_audit_and_clinical_vocabulary_signal_selected_death_modes_2000_2025_soft_clinical_terms_demoted"

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

STAGE_ORDER = [
    "basic_mechanism",
    "preclinical_animal",
    "human_observational",
    "interventional_trial",
    "guideline_ineligible_in_article_corpus",
]

STAGE_LABELS = {
    "basic_mechanism": "Basic\nmechanism",
    "preclinical_animal": "Preclinical\nanimal",
    "human_observational": "Human\nobservational",
    "interventional_trial": "PublicationType-\nsupported trial",
    "guideline_ineligible_in_article_corpus": "Guideline-like\nineligible",
}

STAGE_COLORS = {
    "basic_mechanism": "#DDD1B6",
    "preclinical_animal": "#BFA47D",
    "human_observational": "#8FAFB8",
    "interventional_trial": "#3E6F98",
    "guideline_ineligible_in_article_corpus": "#DADDE2",
}

SHORT_MODE_LABELS = {
    "Ferroptosis": "Ferroptosis",
    "Pyroptosis": "Pyroptosis",
    "NETosis": "NETosis",
    "Necroptosis": "Necroptosis",
    "Immunogenic cell death": "ICD",
    "Cuproptosis": "Cuproptosis",
    "PANoptosis": "PANoptosis",
    "Disulfidptosis": "Disulfidptosis",
}

SOFT_HUMAN_TERMS = ["survival analysis", "serum", "plasma", "biomarker"]
HARD_HUMAN_TERMS = [
    "patient",
    "patients",
    "cohort",
    "diagnostic",
    "prognostic",
    "biopsy",
    "tissue sample",
]


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def ensure_dirs() -> None:
    for subdir in ["pdf", "png", "tiff", "source_data"]:
        (OUT_ROOT / subdir).mkdir(parents=True, exist_ok=True)


def savefig(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_ROOT / "pdf" / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT_ROOT / "png" / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_ROOT / "tiff" / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.5)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def rules_text(df: pd.DataFrame) -> pd.Series:
    cols = [
        "matched_stage_rule",
        "secondary_text_stage_rules",
        "original_matched_stage_rule",
        "audited_pub_type_primary_rules",
    ]
    existing = [c for c in cols if c in df.columns]
    return df[existing].fillna("").agg("|".join, axis=1).str.lower()


def apply_soft_term_reclassification(article: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = article.copy()
    out["audited_stage_original_for_soft_term_variant"] = out["audited_stage"].astype(str)
    rules = rules_text(out)

    soft_pattern = r"human_observational:(?:" + "|".join(re.escape(t) for t in SOFT_HUMAN_TERMS) + r")s?"
    hard_pattern = r"human_observational:(?:" + "|".join(re.escape(t) for t in HARD_HUMAN_TERMS) + r")"

    soft_hit = rules.str.contains(soft_pattern, regex=True, na=False)
    hard_human_hit = rules.str.contains(hard_pattern, regex=True, na=False)
    pubtype_human = (
        rules.str.contains(r"(?:publication_type|original_publication_type):human_observational", regex=True, na=False)
        | out.get("audited_pub_type_primary_stage", pd.Series("", index=out.index))
        .fillna("")
        .astype(str)
        .str.lower()
        .eq("human_observational")
    )
    publication_type_trial = out.get(
        "interventional_supported_by_publication_type_yes_no",
        pd.Series("", index=out.index),
    ).fillna("").astype(str).str.lower().eq("yes")
    preclinical_support = rules.str.contains(r"preclinical_animal:", regex=True, na=False)

    reclassify = (
        out["audited_stage"].astype(str).eq("human_observational")
        & soft_hit
        & ~hard_human_hit
        & ~pubtype_human
        & ~publication_type_trial
        & preclinical_support
    )

    matched_terms = []
    for txt in rules[reclassify]:
        found = [term for term in SOFT_HUMAN_TERMS if re.search(rf"human_observational:{re.escape(term)}s?", txt)]
        matched_terms.append("|".join(found))

    out.loc[reclassify, "audited_stage"] = "preclinical_animal"
    out.loc[reclassify, "audited_tmi_stage"] = "preclinical_animal"
    out.loc[reclassify, "audited_stage_score"] = 0.4
    out.loc[reclassify, "audited_tmi_stage_score"] = 0.4
    out.loc[reclassify, "matched_stage_rule"] = (
        out.loc[reclassify, "matched_stage_rule"].astype(str)
        + "|soft_clinical_terms_demoted_to_preclinical_animal"
    )

    audit = out.loc[
        reclassify,
        [
            "pmid",
            "death_mode",
            "disease_term",
            "title",
            "audited_stage_original_for_soft_term_variant",
            "audited_stage",
            "secondary_text_stage_rules",
            "matched_stage_rule",
        ],
    ].copy()
    audit["soft_terms_triggering_demotion"] = matched_terms
    audit["reclassification_rule"] = (
        "human_observational caused by serum/plasma/biomarker/survival-analysis soft terms; "
        "no hard human-observational or publication-type support; explicit preclinical_animal text support present"
    )

    summary = pd.DataFrame(
        [
            {
                "metric": "records_reclassified_human_observational_to_preclinical_animal",
                "value": int(reclassify.sum()),
            },
            {"metric": "unique_pmids_reclassified", "value": int(out.loc[reclassify, "pmid"].nunique())},
            {"metric": "soft_terms", "value": "|".join(SOFT_HUMAN_TERMS)},
            {"metric": "hard_human_terms_that_block_reclassification", "value": "|".join(HARD_HUMAN_TERMS)},
        ]
    )
    return out, audit, summary


def build_tables(article: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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

    composition_count = (
        article.groupby(["death_mode", "audited_stage"], as_index=False)
        .agg(article_stage_records=("pmid", "size"))
    )
    total_by_mode = article.groupby("death_mode")["pmid"].size().rename("mode_stage_records")
    composition = composition_count.merge(total_by_mode, on="death_mode", how="left")
    composition["stage_share"] = composition["article_stage_records"] / composition["mode_stage_records"]
    complete_rows = []
    for mode in DEATH_ORDER:
        for stage in STAGE_ORDER:
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
    return summary, composition, high_stage


def plot_figure(summary: pd.DataFrame, composition: pd.DataFrame, high_stage: pd.DataFrame, reclassified_n: int) -> None:
    total_records = int(summary.loc[summary["category"].eq("PMID-stage records audited"), "article_stage_records"].iloc[0])
    original_guideline = int(summary.loc[summary["category"].eq("Original guideline assignments"), "article_stage_records"].iloc[0])
    valid_guideline = int(summary.loc[summary["category"].eq("Valid guideline evidence retained"), "article_stage_records"].iloc[0])
    trial_supported = int(summary.loc[summary["category"].eq("PublicationType-supported trials"), "article_stage_records"].iloc[0])

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.8, 3.7),
        constrained_layout=False,
        gridspec_kw={"width_ratios": [1.45, 1.0]},
    )
    y = np.arange(len(DEATH_ORDER))
    left = np.zeros(len(DEATH_ORDER))
    for stage in STAGE_ORDER:
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
        axes[0].barh(y, vals, left=left, height=0.78, color=STAGE_COLORS[stage], label=STAGE_LABELS[stage], **kwargs)
        left += vals
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([SHORT_MODE_LABELS.get(x, x) for x in DEATH_ORDER])
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
    axes[1].set_yticklabels([SHORT_MODE_LABELS.get(x, x) for x in DEATH_ORDER])
    axes[1].set_xlim(0, max(0.82, float(hs["original_high_stage_text_signal_share"].max()) * 1.05))
    axes[1].set_xlabel("Share")
    axes[1].set_title("B. High-stage signal after audit", loc="left", weight="bold")
    axes[1].legend(frameon=False, fontsize=5.7, loc="lower right")
    stats_text = (
        f"{total_records:,} PMID-stage records audited\n"
        f"{reclassified_n:,} soft-term records reclassified to preclinical animal\n"
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
        fontsize=5.55,
        color="#4A5668",
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#CBD5E1", "alpha": 0.92},
    )
    clean_axes(axes[1])
    fig.suptitle(
        "Supplementary Figure S6 analogue. Evidence-stage audit with soft clinical terms demoted",
        x=0.01,
        ha="left",
        fontsize=8.8,
        weight="bold",
    )
    fig.text(
        0.01,
        0.03,
        "Variant rule: serum/plasma/biomarker/survival-analysis text alone is not treated as human observational when stronger human evidence is absent.",
        fontsize=5.75,
        color="#6B7C93",
    )
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.20, top=0.82, wspace=0.32)
    savefig(fig, STEM)


def main() -> int:
    ensure_dirs()
    article = pd.read_csv(ARTICLE_STAGE)
    reclassified, audit, reclass_summary = apply_soft_term_reclassification(article)
    summary, composition, high_stage = build_tables(reclassified)

    audit.to_csv(SRC / f"{STEM}_reclassified_records.csv", index=False)
    reclass_summary.to_csv(SRC / f"{STEM}_reclassification_summary.csv", index=False)
    summary.to_csv(SRC / f"{STEM}_audit_summary.csv", index=False)
    composition.to_csv(SRC / f"{STEM}_evidence_stage_composition.csv", index=False)
    high_stage.to_csv(SRC / f"{STEM}_high_stage_signal_after_audit.csv", index=False)

    plot_figure(summary, composition, high_stage, len(audit))
    print(f"Reclassified records: {len(audit)}")
    print(OUT_ROOT / "pdf" / f"{STEM}.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
