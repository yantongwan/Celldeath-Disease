from __future__ import annotations

import duckdb

from src.queries import (
    get_articles_for_death_mode_disease,
    get_death_mode_disease_summary,
    query_drug_disease,
    query_gene_disease,
)


def test_duplicate_article_stage_rows_removed(conn: duckdb.DuckDBPyConnection) -> None:
    count = conn.execute("SELECT COUNT(*) FROM death_mode_disease_articles").fetchone()[0]
    distinct_count = conn.execute(
        "SELECT COUNT(*) FROM (SELECT DISTINCT pmid, death_mode, disease_term FROM death_mode_disease_articles)"
    ).fetchone()[0]
    assert count == distinct_count


def test_ferroptosis_returns_disease_and_articles(conn: duckdb.DuckDBPyConnection) -> None:
    summary = get_death_mode_disease_summary("Ferroptosis", conn=conn)
    assert not summary.empty
    assert "Breast Neoplasms" in set(summary["disease_term"])
    articles = get_articles_for_death_mode_disease("Ferroptosis", "Breast Neoplasms", conn=conn)
    assert len(articles) == 2
    assert all(articles["pubmed_link"].str.startswith("https://pubmed.ncbi.nlm.nih.gov/"))


def test_gene_casp1_covid_query(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease("CASP1", "COVID-19", conn=conn)
    assert not summary.empty
    assert not articles.empty
    assert set(summary["death_mode"]) == {"Pyroptosis"}


def test_gene_query_supports_disease_system_filter(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease("CASP1", filters={"disease_systems": ["Infections"]}, conn=conn)
    assert not summary.empty
    assert not articles.empty
    assert set(summary["disease_system"]) == {"Infections"}
    assert set(articles["disease_system"]) == {"Infections"}


def test_gene_query_supports_multiple_selected_symbols(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease(gene_symbols=["ACSL4", "CASP1", "ASCL4"], conn=conn)
    assert set(summary["gene_symbol"]) == {"ACSL4", "CASP1"}
    assert set(articles["gene_symbol"]) == {"ACSL4", "CASP1"}


def test_non_existing_gene_returns_empty_valid_tables(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease("NO_SUCH_GENE", "COVID-19", conn=conn)
    assert summary.empty
    assert articles.empty
    assert "gene_symbol" in summary.columns
    assert "pmid" in articles.columns


def test_ascl4_hidden_by_default(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_gene_disease("ASCL4", show_artifacts=False, conn=conn)
    assert summary.empty
    assert articles.empty
    visible_summary, visible_articles = query_gene_disease("ASCL4", show_artifacts=True, conn=conn)
    assert not visible_summary.empty
    assert not visible_articles.empty
    assert bool(visible_summary["artifact_flag"].iloc[0]) is True


def test_drug_cisplatin_disease_query(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drug_disease("cisplatin", "Breast", conn=conn)
    assert not summary.empty
    assert not articles.empty
    assert set(summary["death_mode"]) == {"Ferroptosis"}


def test_drug_query_supports_multiple_selected_names(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drug_disease(drug_names=["cisplatin", "missing drug"], disease_search="Breast", conn=conn)
    assert not summary.empty
    assert not articles.empty
    assert set(summary["normalized_drug_name"]) == {"cisplatin"}
    assert set(articles["normalized_drug_name"]) == {"cisplatin"}


def test_drug_query_supports_disease_system_filter(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drug_disease("cisplatin", disease_systems=["Neoplasms"], conn=conn)
    assert not summary.empty
    assert not articles.empty
    assert set(summary["disease_system"]) == {"Neoplasms"}
    assert set(articles["disease_system"]) == {"Neoplasms"}


def test_negated_drug_relations_excluded_by_default(conn: duckdb.DuckDBPyConnection) -> None:
    summary, articles = query_drug_disease("cisplatin", "Lung", exclude_negated=True, conn=conn)
    assert summary.empty
    assert articles.empty
    visible_summary, visible_articles = query_drug_disease(
        "cisplatin", "Lung", exclude_negated=False, conn=conn
    )
    assert not visible_summary.empty
    assert not visible_articles.empty
