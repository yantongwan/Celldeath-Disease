from __future__ import annotations

from collections import OrderedDict
from typing import Iterable

import duckdb
import pandas as pd


TABLE_SCHEMAS: dict[str, OrderedDict[str, str]] = {
    "articles": OrderedDict(
        [
            ("pmid", "VARCHAR PRIMARY KEY"),
            ("title", "VARCHAR"),
            ("journal", "VARCHAR"),
            ("publication_year", "INTEGER"),
            ("publication_date", "VARCHAR"),
            ("doi", "VARCHAR"),
            ("abstract", "VARCHAR"),
            ("publication_types", "VARCHAR"),
            ("country", "VARCHAR"),
            ("evidence_stage", "VARCHAR"),
            ("high_stage_text_signal", "BOOLEAN"),
            ("publicationtype_supported_trial", "BOOLEAN"),
            ("is_oncology", "BOOLEAN"),
            ("tumor_family", "VARCHAR"),
        ]
    ),
    "death_mode_disease_articles": OrderedDict(
        [
            ("pmid", "VARCHAR"),
            ("death_mode", "VARCHAR"),
            ("disease_term", "VARCHAR"),
            ("disease_system", "VARCHAR"),
            ("is_neoplasm", "BOOLEAN"),
            ("tumor_family", "VARCHAR"),
            ("evidence_stage", "VARCHAR"),
        ]
    ),
    "death_mode_disease_pairs": OrderedDict(
        [
            ("death_mode", "VARCHAR"),
            ("disease_term", "VARCHAR"),
            ("disease_system", "VARCHAR"),
            ("is_neoplasm", "BOOLEAN"),
            ("tumor_family", "VARCHAR"),
            ("pair_count", "INTEGER"),
            ("unique_pmids", "INTEGER"),
            ("first_year", "INTEGER"),
            ("latest_year", "INTEGER"),
            ("literature_hub_score", "DOUBLE"),
            ("selectivity_score", "DOUBLE"),
        ]
    ),
    "gene_mentions": OrderedDict(
        [
            ("pmid", "VARCHAR"),
            ("gene_symbol", "VARCHAR"),
            ("gene_alias", "VARCHAR"),
            ("source_field", "VARCHAR"),
            ("matched_text", "VARCHAR"),
            ("context_sentence", "VARCHAR"),
            ("artifact_flag", "BOOLEAN"),
        ]
    ),
    "gene_death_disease_articles": OrderedDict(
        [
            ("pmid", "VARCHAR"),
            ("gene_symbol", "VARCHAR"),
            ("death_mode", "VARCHAR"),
            ("disease_term", "VARCHAR"),
            ("disease_system", "VARCHAR"),
            ("is_oncology", "BOOLEAN"),
            ("tumor_family", "VARCHAR"),
            ("source_field", "VARCHAR"),
            ("context_sentence", "VARCHAR"),
            ("evidence_stage", "VARCHAR"),
            ("artifact_flag", "BOOLEAN"),
        ]
    ),
    "drug_mentions": OrderedDict(
        [
            ("pmid", "VARCHAR"),
            ("drug_name", "VARCHAR"),
            ("normalized_drug_name", "VARCHAR"),
            ("reference_class_label", "VARCHAR"),
            ("approved_clinical_drug_flag", "BOOLEAN"),
            ("fda_approved_flag", "BOOLEAN"),
            ("source", "VARCHAR"),
            ("relation_cue", "VARCHAR"),
            ("relation_sentence", "VARCHAR"),
            ("negation_flag", "BOOLEAN"),
            ("token_distance", "INTEGER"),
            ("relation_confidence", "VARCHAR"),
        ]
    ),
    "drug_death_disease_articles": OrderedDict(
        [
            ("pmid", "VARCHAR"),
            ("normalized_drug_name", "VARCHAR"),
            ("death_mode", "VARCHAR"),
            ("disease_term", "VARCHAR"),
            ("tumor_family", "VARCHAR"),
            ("reference_class_label", "VARCHAR"),
            ("approved_clinical_drug_flag", "BOOLEAN"),
            ("fda_approved_flag", "BOOLEAN"),
            ("relation_cue", "VARCHAR"),
            ("relation_sentence", "VARCHAR"),
            ("negation_flag", "BOOLEAN"),
            ("token_distance", "INTEGER"),
            ("evidence_stage", "VARCHAR"),
        ]
    ),
    "synonyms": OrderedDict(
        [
            ("entity_type", "VARCHAR"),
            ("canonical_name", "VARCHAR"),
            ("synonym", "VARCHAR"),
        ]
    ),
    "build_sources": OrderedDict(
        [
            ("logical_table", "VARCHAR"),
            ("source_file", "VARCHAR"),
            ("rows_loaded", "INTEGER"),
            ("column_mapping", "VARCHAR"),
        ]
    ),
    "validation_warnings": OrderedDict(
        [
            ("severity", "VARCHAR"),
            ("message", "VARCHAR"),
        ]
    ),
}


VIEW_NAMES = (
    "v_death_mode_summary",
    "v_disease_by_death_mode",
    "v_gene_query_summary",
    "v_drug_query_summary",
)


def empty_frame(table: str) -> pd.DataFrame:
    return pd.DataFrame(columns=list(TABLE_SCHEMAS[table].keys()))


def ensure_columns(df: pd.DataFrame, table: str) -> pd.DataFrame:
    schema = TABLE_SCHEMAS[table]
    out = df.copy()
    for col in schema:
        if col not in out.columns:
            out[col] = None
    return out[list(schema.keys())]


def _column_definitions(schema: OrderedDict[str, str]) -> str:
    return ", ".join(f"{name} {sql_type}" for name, sql_type in schema.items())


def create_tables(conn: duckdb.DuckDBPyConnection, drop: bool = False) -> None:
    if drop:
        for view in VIEW_NAMES:
            conn.execute(f"DROP VIEW IF EXISTS {view}")
        for table in TABLE_SCHEMAS:
            conn.execute(f"DROP TABLE IF EXISTS {table}")
    for table, schema in TABLE_SCHEMAS.items():
        conn.execute(f"CREATE TABLE IF NOT EXISTS {table} ({_column_definitions(schema)})")


def insert_dataframe(conn: duckdb.DuckDBPyConnection, table: str, df: pd.DataFrame) -> None:
    if df.empty:
        return
    prepared = ensure_columns(df, table)
    conn.register("_insert_df", prepared)
    cols = ", ".join(prepared.columns)
    conn.execute(f"INSERT INTO {table} ({cols}) SELECT {cols} FROM _insert_df")
    conn.unregister("_insert_df")


def create_views(conn: duckdb.DuckDBPyConnection) -> None:
    conn.execute(
        """
        CREATE OR REPLACE VIEW v_death_mode_summary AS
        WITH pair_counts AS (
            SELECT
                death_mode,
                disease_term,
                COUNT(DISTINCT pmid) AS unique_pmids
            FROM death_mode_disease_articles
            GROUP BY death_mode, disease_term
        ),
        ranked AS (
            SELECT
                death_mode,
                disease_term,
                unique_pmids,
                ROW_NUMBER() OVER (
                    PARTITION BY death_mode
                    ORDER BY unique_pmids DESC, disease_term
                ) AS rn
            FROM pair_counts
        ),
        top_diseases AS (
            SELECT
                death_mode,
                STRING_AGG(
                    disease_term || ' (' || CAST(unique_pmids AS VARCHAR) || ')',
                    '; ' ORDER BY unique_pmids DESC, disease_term
                ) AS top_diseases
            FROM ranked
            WHERE rn <= 5
            GROUP BY death_mode
        ),
        approved_drugs AS (
            SELECT
                death_mode,
                COUNT(DISTINCT normalized_drug_name) AS approved_drug_count
            FROM drug_death_disease_articles
            WHERE approved_clinical_drug_flag IS TRUE OR fda_approved_flag IS TRUE
            GROUP BY death_mode
        )
        SELECT
            d.death_mode,
            COUNT(DISTINCT d.pmid) AS total_unique_pmids,
            COUNT(DISTINCT d.disease_term) AS disease_breadth,
            COALESCE(t.top_diseases, '') AS top_diseases,
            COUNT(DISTINCT CASE
                WHEN a.publicationtype_supported_trial IS TRUE THEN d.pmid
                ELSE NULL
            END) AS trial_supported_count,
            AVG(CASE WHEN COALESCE(a.is_oncology, d.is_neoplasm) IS TRUE THEN 1.0 ELSE 0.0 END)
                AS oncology_fraction,
            COALESCE(ad.approved_drug_count, 0) AS approved_drug_count
        FROM death_mode_disease_articles d
        LEFT JOIN articles a ON a.pmid = d.pmid
        LEFT JOIN top_diseases t ON t.death_mode = d.death_mode
        LEFT JOIN approved_drugs ad ON ad.death_mode = d.death_mode
        GROUP BY d.death_mode, t.top_diseases, ad.approved_drug_count
        """
    )
    conn.execute(
        """
        CREATE OR REPLACE VIEW v_disease_by_death_mode AS
        SELECT
            death_mode,
            disease_term,
            pair_count,
            unique_pmids,
            first_year,
            latest_year,
            disease_system,
            literature_hub_score,
            selectivity_score
        FROM death_mode_disease_pairs
        """
    )
    conn.execute(
        """
        CREATE OR REPLACE VIEW v_gene_query_summary AS
        WITH stage_counts AS (
            SELECT
                gene_symbol,
                disease_term,
                death_mode,
                COALESCE(evidence_stage, 'missing') AS evidence_stage,
                COUNT(DISTINCT pmid) AS stage_pmids
            FROM gene_death_disease_articles
            WHERE artifact_flag IS NOT TRUE
            GROUP BY gene_symbol, disease_term, death_mode, COALESCE(evidence_stage, 'missing')
        ),
        stage_text AS (
            SELECT
                gene_symbol,
                disease_term,
                death_mode,
                STRING_AGG(
                    evidence_stage || '=' || CAST(stage_pmids AS VARCHAR),
                    '; ' ORDER BY evidence_stage
                ) AS evidence_stage_distribution
            FROM stage_counts
            GROUP BY gene_symbol, disease_term, death_mode
        )
        SELECT
            g.gene_symbol,
            g.disease_term,
            g.death_mode,
            COUNT(DISTINCT g.pmid) AS unique_pmids,
            MIN(a.publication_year) AS first_year,
            MAX(a.publication_year) AS latest_year,
            COALESCE(s.evidence_stage_distribution, '') AS evidence_stage_distribution,
            AVG(CASE WHEN g.is_oncology IS TRUE THEN 1.0 ELSE 0.0 END) AS oncology_fraction
        FROM gene_death_disease_articles g
        LEFT JOIN articles a ON a.pmid = g.pmid
        LEFT JOIN stage_text s
            ON s.gene_symbol = g.gene_symbol
            AND s.disease_term = g.disease_term
            AND s.death_mode = g.death_mode
        WHERE g.artifact_flag IS NOT TRUE
        GROUP BY g.gene_symbol, g.disease_term, g.death_mode, s.evidence_stage_distribution
        """
    )
    conn.execute(
        """
        CREATE OR REPLACE VIEW v_drug_query_summary AS
        SELECT
            d.normalized_drug_name,
            d.disease_term,
            d.death_mode,
            COUNT(DISTINCT d.pmid) AS unique_pmids,
            BOOL_OR(d.approved_clinical_drug_flag) AS approved_clinical_drug_flag,
            STRING_AGG(DISTINCT COALESCE(d.reference_class_label, ''), '; ') AS reference_class_label,
            MIN(a.publication_year) AS first_year,
            MAX(a.publication_year) AS latest_year,
            STRING_AGG(DISTINCT COALESCE(d.relation_cue, ''), '; ') AS relation_cues
        FROM drug_death_disease_articles d
        LEFT JOIN articles a ON a.pmid = d.pmid
        WHERE d.negation_flag IS NOT TRUE
        GROUP BY d.normalized_drug_name, d.disease_term, d.death_mode
        """
    )


def table_names(include_internal: bool = False) -> Iterable[str]:
    for table in TABLE_SCHEMAS:
        if include_internal or table not in {"build_sources", "validation_warnings"}:
            yield table
