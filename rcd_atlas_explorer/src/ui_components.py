from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st

from .io_utils import load_yaml


def load_ui_config(project_root: Path) -> dict:
    return load_yaml(project_root / "config" / "ui_config.yaml")


def apply_page_style() -> None:
    st.markdown(
        """
        <style>
        .block-container {padding-top: 1.4rem; padding-bottom: 2.5rem;}
        div[data-testid="stMetric"] {
            background: #f8fafc;
            border: 1px solid #e5e7eb;
            border-radius: 8px;
            padding: 0.7rem 0.8rem;
        }
        div[data-testid="stAlert"] {border-radius: 8px;}
        .small-note {color: #475569; font-size: 0.9rem;}
        </style>
        """,
        unsafe_allow_html=True,
    )


def csv_bytes(df: pd.DataFrame) -> bytes:
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue().encode("utf-8")


def download_table(df: pd.DataFrame, label: str, filename: str, key: str) -> None:
    st.download_button(
        label,
        data=csv_bytes(df),
        file_name=filename,
        mime="text/csv",
        key=key,
        disabled=df.empty,
    )


def pmid_text_area(df: pd.DataFrame, key: str) -> None:
    pmids = []
    if "pmid" in df:
        pmids = sorted({str(value) for value in df["pmid"].dropna() if str(value).strip()})
    st.text_area("PMID list", "\n".join(pmids), height=110, key=key)


def show_dataframe(df: pd.DataFrame, key: str, height: int = 420) -> None:
    column_config = {}
    for col in ["pubmed_link", "doi_link"]:
        if col in df.columns:
            column_config[col] = st.column_config.LinkColumn(col)
    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        height=height,
        column_config=column_config,
        key=key,
    )


def safe_bar_chart(data: pd.Series | pd.DataFrame, key: str) -> None:
    try:
        st.bar_chart(data)
    except Exception as exc:
        st.warning(
            "Chart rendering is unavailable in this Python environment. "
            f"Showing the same values as a table instead. Error: {exc}"
        )
        if isinstance(data, pd.Series):
            fallback = data.rename("value").reset_index()
        else:
            fallback = data.reset_index()
        st.dataframe(fallback, use_container_width=True, hide_index=True, key=key)


def metric_row(metrics: dict[str, int]) -> None:
    cols = st.columns(len(metrics))
    for col, (label, value) in zip(cols, metrics.items()):
        col.metric(label, f"{value:,}")


def summary_cards_from_articles(df: pd.DataFrame) -> None:
    if df.empty:
        metric_row(
            {
                "Articles": 0,
                "Death modes": 0,
                "Diseases": 0,
                "First year": 0,
                "Latest year": 0,
            }
        )
        return
    years = pd.to_numeric(df.get("publication_year", pd.Series(dtype=float)), errors="coerce")
    metric_row(
        {
            "Articles": int(df["pmid"].nunique()) if "pmid" in df else len(df),
            "Death modes": int(df["death_mode"].nunique()) if "death_mode" in df else 0,
            "Diseases": int(df["disease_term"].nunique()) if "disease_term" in df else 0,
            "First year": int(years.min()) if years.notna().any() else 0,
            "Latest year": int(years.max()) if years.notna().any() else 0,
        }
    )


DATA_DICTIONARY = """
# Data Dictionary

## Core article columns

- `pmid`: PubMed identifier stored as a string.
- `title`, `journal`, `publication_year`, `publication_date`, `doi`, `abstract`: article metadata.
- `publication_types`: PubMed publication-type text when available.
- `evidence_stage`: source-level stage label. This is descriptive and not a validation grade.
- `publicationtype_supported_trial`: true when publication-type text supports trial-level evidence.
- `is_oncology`, `tumor_family`: oncology-oriented flags from source data or disease-term inference.

## Atlas association columns

- `death_mode`: one of the eight canonical regulated cell-death modes.
- `disease_term`: normalized disease term from source tables.
- `disease_system`: optional disease-system grouping.
- `pair_count`, `unique_pmids`: literature-volume measures for a death-mode disease pair.
- `literature_hub_score`, `selectivity_score`: manuscript-derived descriptive metrics when available.

## Gene columns

- `gene_symbol`: uppercase gene symbol after normalization.
- `source_field`, `matched_text`, `context_sentence`: text-mining provenance fields.
- `artifact_flag`: true for entries excluded from biological summaries by default. `ASCL4` in ferroptosis contexts is flagged as a likely `ACSL4` spelling/entity-normalization artefact.

## Drug columns

- `normalized_drug_name`: conservative normalized drug label.
- `reference_class_label`: source drug/reference-class label where available.
- `approved_clinical_drug_flag`, `fda_approved_flag`: source flags; absence does not imply non-approval.
- `relation_cue`, `relation_sentence`, `negation_flag`, `token_distance`: sentence-level relation provenance.

## Interpretation boundary

This application summarizes literature-level associations from PubMed-derived data. A death-mode-disease pair, gene co-mention or drug-death relation candidate does not prove biological causality, pathway activation, therapeutic efficacy or clinical validity.
"""
