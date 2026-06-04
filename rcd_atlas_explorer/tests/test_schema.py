from __future__ import annotations

import duckdb

from src.schema import TABLE_SCHEMAS


def test_required_tables_exist(conn: duckdb.DuckDBPyConnection) -> None:
    tables = {
        row[0]
        for row in conn.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
        ).fetchall()
    }
    for table in TABLE_SCHEMAS:
        assert table in tables


def test_pmids_are_stored_as_text(conn: duckdb.DuckDBPyConnection) -> None:
    pmids = conn.execute("SELECT pmid FROM articles ORDER BY pmid").fetchall()
    assert pmids
    assert all(isinstance(row[0], str) for row in pmids)
