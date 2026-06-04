from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd


def connect(db_path: str | Path, read_only: bool = True) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(str(db_path), read_only=read_only)


def table_count(conn: duckdb.DuckDBPyConnection, table: str) -> int:
    return int(conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])


def table_counts(conn: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    rows = []
    for table in [
        "articles",
        "death_mode_disease_articles",
        "death_mode_disease_pairs",
        "gene_mentions",
        "gene_death_disease_articles",
        "drug_mentions",
        "drug_death_disease_articles",
        "synonyms",
    ]:
        rows.append({"table": table, "rows": table_count(conn, table)})
    return pd.DataFrame(rows)


def overview_metrics(conn: duckdb.DuckDBPyConnection) -> dict[str, int]:
    queries = {
        "unique_pmids": "SELECT COUNT(DISTINCT pmid) FROM articles",
        "death_modes": "SELECT COUNT(DISTINCT death_mode) FROM death_mode_disease_articles",
        "diseases": "SELECT COUNT(DISTINCT disease_term) FROM death_mode_disease_articles",
        "death_mode_disease_pairs": (
            "SELECT COUNT(*) FROM (SELECT DISTINCT death_mode, disease_term "
            "FROM death_mode_disease_articles)"
        ),
        "genes": "SELECT COUNT(DISTINCT gene_symbol) FROM gene_mentions",
        "drugs": "SELECT COUNT(DISTINCT normalized_drug_name) FROM drug_mentions",
    }
    return {name: int(conn.execute(sql).fetchone()[0] or 0) for name, sql in queries.items()}


def year_bounds(conn: duckdb.DuckDBPyConnection) -> tuple[int, int]:
    row = conn.execute(
        "SELECT MIN(publication_year), MAX(publication_year) FROM articles WHERE publication_year IS NOT NULL"
    ).fetchone()
    start = int(row[0] or 2000)
    end = int(row[1] or start)
    return start, end


def distinct_values(conn: duckdb.DuckDBPyConnection, table: str, column: str) -> list[str]:
    rows = conn.execute(
        f"SELECT DISTINCT {column} FROM {table} WHERE {column} IS NOT NULL AND {column} <> '' ORDER BY {column}"
    ).fetchall()
    return [str(row[0]) for row in rows]
