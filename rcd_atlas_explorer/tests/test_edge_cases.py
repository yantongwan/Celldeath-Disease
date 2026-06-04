from __future__ import annotations

import duckdb

from src.queries import (
    get_articles_for_death_mode_disease,
    global_article_search,
    query_drug_disease,
    query_drugs_for_death_mode_disease,
    query_gene_disease,
)


def test_empty_death_or_disease_inputs_return_valid_tables(conn: duckdb.DuckDBPyConnection) -> None:
    summary = get_articles_for_death_mode_disease("", "Breast Neoplasms", conn=conn)
    assert summary.empty
    assert {"pmid", "pubmed_link", "abstract_snippet"}.issubset(summary.columns)


def test_relation_cue_filter_limits_drug_results(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drug_disease(
        "cisplatin",
        "Breast",
        relation_cues=["induces"],
        conn=conn,
    )
    assert not summary.empty
    assert not articles.empty
    assert set(articles["relation_cue"]) == {"induces"}

    missing_summary, missing_articles = query_drug_disease(
        "cisplatin",
        "Breast",
        relation_cues=["not associated"],
        conn=conn,
    )
    assert missing_summary.empty
    assert missing_articles.empty


def test_pair_drug_query_supports_relation_cue_filter(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drugs_for_death_mode_disease(
        "Ferroptosis",
        "Breast Neoplasms",
        relation_cues=["induces"],
        conn=conn,
    )
    assert not summary.empty
    assert not articles.empty
    assert set(summary["normalized_drug_name"]) == {"cisplatin"}


def test_global_search_without_text_can_filter_to_gene_or_drug(conn: duckdb.DuckDBPyConnection) -> None:
    gene_rows = global_article_search(filters={"gene_present": True}, conn=conn)
    drug_rows = global_article_search(filters={"drug_present": True}, conn=conn)
    assert not gene_rows.empty
    assert not drug_rows.empty
    assert "genes" in gene_rows.columns
    assert "drugs" in drug_rows.columns


def test_gene_query_with_empty_inputs_returns_valid_tables(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease(conn=conn)
    assert summary.empty
    assert articles.empty
    assert {"gene_symbol", "artifact_flag"}.issubset(summary.columns)
    assert {"pmid", "pubmed_link"}.issubset(articles.columns)
