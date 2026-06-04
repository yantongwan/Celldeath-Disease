from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .normalization import CANONICAL_DEATH_MODES, normalize_pmid
from .schema import TABLE_SCHEMAS


def validate_tables(
    tables: dict[str, pd.DataFrame],
    source_records: list[dict[str, Any]],
    warnings: list[str],
    year_min: int = 2000,
    year_max: int = 2025,
) -> tuple[list[dict[str, str]], str]:
    records: list[dict[str, str]] = [{"severity": "warning", "message": msg} for msg in warnings]

    for table, schema in TABLE_SCHEMAS.items():
        if table in {"build_sources", "validation_warnings"}:
            continue
        df = tables.get(table, pd.DataFrame())
        missing_cols = [col for col in schema if col not in df.columns]
        if missing_cols:
            records.append(
                {
                    "severity": "error",
                    "message": f"{table}: missing final columns {missing_cols}",
                }
            )
        if "pmid" in df.columns:
            bad_pmids = df["pmid"].map(normalize_pmid).eq("").sum()
            if bad_pmids:
                records.append(
                    {
                        "severity": "warning",
                        "message": f"{table}: {bad_pmids} rows have missing or unparseable PMID",
                    }
                )
        if "death_mode" in df.columns and not df.empty:
            unmapped = sorted(
                {
                    str(value)
                    for value in df["death_mode"].dropna().unique()
                    if value and value not in CANONICAL_DEATH_MODES
                }
            )
            if unmapped:
                records.append(
                    {
                        "severity": "warning",
                        "message": f"{table}: unmapped death modes {unmapped[:20]}",
                    }
                )
        if "publication_year" in df.columns and not df.empty:
            years = pd.to_numeric(df["publication_year"], errors="coerce")
            outside = years.notna() & ((years < year_min) | (years > year_max))
            if int(outside.sum()):
                records.append(
                    {
                        "severity": "warning",
                        "message": (
                            f"{table}: {int(outside.sum())} publication years outside "
                            f"{year_min}-{year_max}"
                        ),
                    }
                )

    articles = tables.get("articles", pd.DataFrame())
    if not articles.empty:
        for col in ["title", "journal", "publication_year"]:
            if col in articles:
                missing = articles[col].isna() | articles[col].astype(str).str.strip().eq("")
                records.append(
                    {
                        "severity": "info",
                        "message": f"articles: {int(missing.sum())} rows missing {col}",
                    }
                )

    dmd = tables.get("death_mode_disease_articles", pd.DataFrame())
    if not dmd.empty:
        duplicate_cols = ["pmid", "death_mode", "disease_term"]
        duplicates = dmd.duplicated(subset=duplicate_cols).sum()
        if int(duplicates):
            records.append(
                {
                    "severity": "warning",
                    "message": (
                        "death_mode_disease_articles: "
                        f"{int(duplicates)} duplicate PMID-death-disease rows"
                    ),
                }
            )

    genes = tables.get("gene_death_disease_articles", pd.DataFrame())
    if not genes.empty and "artifact_flag" in genes:
        flagged = genes["artifact_flag"].fillna(False).astype(bool).sum()
        records.append(
            {
                "severity": "info",
                "message": f"gene_death_disease_articles: {int(flagged)} ASCL4/artifact rows flagged",
            }
        )

    drugs = tables.get("drug_death_disease_articles", pd.DataFrame())
    if not drugs.empty and "negation_flag" in drugs:
        negated = drugs["negation_flag"].fillna(False).astype(bool).sum()
        records.append(
            {
                "severity": "info",
                "message": f"drug_death_disease_articles: {int(negated)} negated rows retained",
            }
        )

    report = render_validation_report(tables, source_records, records)
    return records, report


def render_validation_report(
    tables: dict[str, pd.DataFrame],
    source_records: list[dict[str, Any]],
    validation_records: list[dict[str, str]],
) -> str:
    lines: list[str] = [
        "# Validation Report",
        "",
        "## Source Files Loaded",
        "",
    ]
    if source_records:
        lines.extend(["| logical_table | rows | source_file | column_mapping |", "|---|---:|---|---|"])
        for record in source_records:
            lines.append(
                f"| {record['logical_table']} | {record['rows_loaded']} | "
                f"{Path(record['source_file']).name} | `{record['column_mapping']}` |"
            )
    else:
        lines.append("No raw source files were loaded.")

    lines.extend(["", "## Final Table Row Counts", ""])
    lines.extend(["| table | rows | columns |", "|---|---:|---:|"])
    for table in sorted(tables):
        df = tables[table]
        lines.append(f"| {table} | {len(df)} | {len(df.columns)} |")

    lines.extend(["", "## Warnings and Notes", ""])
    if validation_records:
        lines.extend(["| severity | message |", "|---|---|"])
        for record in validation_records:
            message = str(record["message"]).replace("|", "\\|")
            lines.append(f"| {record['severity']} | {message} |")
    else:
        lines.append("No validation warnings were generated.")

    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "Drug-death-disease associations, gene co-mentions and disease-pair counts are "
            "literature-level signals. They do not establish biological causality, pathway "
            "activation, therapeutic efficacy or clinical validity.",
            "",
        ]
    )
    return "\n".join(lines)


def write_validation_report(path: Path, report: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report, encoding="utf-8")
