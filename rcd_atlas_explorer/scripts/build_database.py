from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.io_utils import load_logical_table, load_yaml
from src.normalization import (
    CANONICAL_DEATH_MODES,
    as_bool,
    as_float,
    as_int,
    flag_gene_artifact,
    infer_oncology_flag,
    normalize_death_mode,
    normalize_disease,
    normalize_drug_name,
    normalize_gene_symbol,
    normalize_pmid,
)
from src.schema import (
    TABLE_SCHEMAS,
    create_tables,
    create_views,
    empty_frame,
    ensure_columns,
    insert_dataframe,
    table_names,
)
from src.validation import validate_tables, write_validation_report


BOOL_COLUMNS = {
    "articles": ["high_stage_text_signal", "publicationtype_supported_trial", "is_oncology"],
    "death_mode_disease_articles": ["is_neoplasm"],
    "death_mode_disease_pairs": ["is_neoplasm"],
    "gene_mentions": ["artifact_flag"],
    "gene_death_disease_articles": ["is_oncology", "artifact_flag"],
    "drug_mentions": ["approved_clinical_drug_flag", "fda_approved_flag", "negation_flag"],
    "drug_death_disease_articles": [
        "approved_clinical_drug_flag",
        "fda_approved_flag",
        "negation_flag",
    ],
}

INT_COLUMNS = {
    "articles": ["publication_year"],
    "death_mode_disease_pairs": ["pair_count", "unique_pmids", "first_year", "latest_year"],
    "drug_mentions": ["token_distance"],
    "drug_death_disease_articles": ["token_distance"],
}

FLOAT_COLUMNS = {
    "death_mode_disease_pairs": ["literature_hub_score", "selectivity_score"],
}


def _first_nonempty(series: pd.Series) -> Any:
    for value in series:
        if pd.notna(value) and str(value).strip() != "":
            return value
    return None


def _collapse_by_key(df: pd.DataFrame, key: str, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    grouped = df.groupby(key, dropna=False)[columns].agg(_first_nonempty).reset_index()
    return grouped


def _bool_series(series: pd.Series) -> pd.Series:
    return series.map(as_bool).astype("boolean")


def _int_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.map(as_int), errors="coerce").astype("Int64")


def _float_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.map(as_float), errors="coerce")


def _is_approved_drug_label(value: Any) -> bool:
    text = str(value or "").strip().lower().replace("_", " ")
    return "approved clinical drug" in text or "fda approved" in text


def _infer_drug_approval_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    label_flag = out.get("reference_class_label", pd.Series("", index=out.index)).map(
        _is_approved_drug_label
    )
    for col in ["approved_clinical_drug_flag", "fda_approved_flag"]:
        explicit = (
            out[col].map(as_bool)
            if col in out
            else pd.Series([None] * len(out), index=out.index)
        )
        out[col] = [value if value is not None else inferred for value, inferred in zip(explicit, label_flag)]
    return out


def _prefer_approved_rows(df: pd.DataFrame) -> pd.DataFrame:
    sort_cols = [
        col
        for col in ["approved_clinical_drug_flag", "fda_approved_flag"]
        if col in df.columns
    ]
    if not sort_cols:
        return df
    return df.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last")


def cast_for_schema(df: pd.DataFrame, table: str) -> pd.DataFrame:
    out = df.copy()
    for col in TABLE_SCHEMAS[table]:
        if col not in out:
            out[col] = None
    for col in BOOL_COLUMNS.get(table, []):
        out[col] = _bool_series(out[col])
    for col in INT_COLUMNS.get(table, []):
        out[col] = _int_series(out[col])
    for col in FLOAT_COLUMNS.get(table, []):
        out[col] = _float_series(out[col])
    text_cols = [
        col
        for col, sql_type in TABLE_SCHEMAS[table].items()
        if sql_type.startswith("VARCHAR")
    ]
    for col in text_cols:
        out[col] = out[col].fillna("").astype(str)
    return ensure_columns(out, table)


def normalize_articles(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("articles")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    out = out[out["pmid"] != ""].copy()
    if "publication_year" in out:
        out["publication_year"] = out["publication_year"].map(as_int)
    else:
        out["publication_year"] = None
    if "publicationtype_supported_trial" not in out:
        out["publicationtype_supported_trial"] = None
    publication_types = out.get("publication_types", pd.Series("", index=out.index)).astype(str)
    inferred_trial = publication_types.str.contains("clinical trial|randomized|trial", case=False, regex=True)
    out["publicationtype_supported_trial"] = out["publicationtype_supported_trial"].map(as_bool)
    out.loc[out["publicationtype_supported_trial"].isna(), "publicationtype_supported_trial"] = inferred_trial
    if "is_oncology" in out:
        out["is_oncology"] = out["is_oncology"].map(as_bool)
    else:
        out["is_oncology"] = None
    if "high_stage_text_signal" in out:
        out["high_stage_text_signal"] = out["high_stage_text_signal"].map(as_bool)
    else:
        out["high_stage_text_signal"] = None
    for col in TABLE_SCHEMAS["articles"]:
        if col not in out:
            out[col] = None
    collapsed = _collapse_by_key(out, "pmid", [col for col in TABLE_SCHEMAS["articles"] if col != "pmid"])
    return cast_for_schema(collapsed, "articles")


def normalize_death_mode_disease_articles(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("death_mode_disease_articles")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    out["death_mode"] = out.get("death_mode", "").map(normalize_death_mode)
    out["disease_term"] = out.get("disease_term", "").map(normalize_disease)
    out = out[(out["pmid"] != "") & (out["death_mode"] != "") & (out["disease_term"] != "")].copy()
    out = out[out["death_mode"].isin(CANONICAL_DEATH_MODES)].copy()
    if "is_neoplasm" in out:
        out["is_neoplasm"] = [
            infer_oncology_flag(disease, explicit)
            for disease, explicit in zip(out["disease_term"], out["is_neoplasm"])
        ]
    else:
        out["is_neoplasm"] = out["disease_term"].map(infer_oncology_flag)
    for col in TABLE_SCHEMAS["death_mode_disease_articles"]:
        if col not in out:
            out[col] = None
    out = out.drop_duplicates(subset=["pmid", "death_mode", "disease_term"])
    return cast_for_schema(out, "death_mode_disease_articles")


def normalize_pairs(df: pd.DataFrame, article_stage: pd.DataFrame, articles: pd.DataFrame) -> pd.DataFrame:
    if df.empty and article_stage.empty:
        return empty_frame("death_mode_disease_pairs")
    if df.empty:
        base = article_stage.groupby(["death_mode", "disease_term"], dropna=False).agg(
            unique_pmids=("pmid", "nunique"),
            pair_count=("pmid", "nunique"),
            disease_system=("disease_system", _first_nonempty),
            is_neoplasm=("is_neoplasm", _first_nonempty),
            tumor_family=("tumor_family", _first_nonempty),
        )
        out = base.reset_index()
    else:
        out = df.copy()
        out["death_mode"] = out.get("death_mode", "").map(normalize_death_mode)
        out["disease_term"] = out.get("disease_term", "").map(normalize_disease)
        out = out[(out["death_mode"] != "") & (out["disease_term"] != "")].copy()
        out = out[out["death_mode"].isin(CANONICAL_DEATH_MODES)].copy()
    if not article_stage.empty:
        year_stage = article_stage[["pmid", "death_mode", "disease_term"]].merge(
            articles[["pmid", "publication_year"]], on="pmid", how="left"
        )
        year_summary = (
            year_stage.groupby(["death_mode", "disease_term"], dropna=False)
            .agg(
                computed_unique_pmids=("pmid", "nunique"),
                computed_first_year=("publication_year", "min"),
                computed_latest_year=("publication_year", "max"),
            )
            .reset_index()
        )
        out = out.merge(year_summary, on=["death_mode", "disease_term"], how="left")
    else:
        out["computed_unique_pmids"] = None
        out["computed_first_year"] = None
        out["computed_latest_year"] = None
    for col in ["pair_count", "unique_pmids"]:
        if col not in out:
            out[col] = None
    out["unique_pmids"] = out["unique_pmids"].replace("", pd.NA).fillna(out["computed_unique_pmids"])
    out["pair_count"] = out["pair_count"].replace("", pd.NA).fillna(out["unique_pmids"])
    if "first_year" not in out:
        out["first_year"] = None
    if "latest_year" not in out:
        out["latest_year"] = None
    out["first_year"] = out["first_year"].replace("", pd.NA).fillna(out["computed_first_year"])
    out["latest_year"] = out["latest_year"].replace("", pd.NA).fillna(out["computed_latest_year"])
    if "is_neoplasm" in out:
        out["is_neoplasm"] = [
            infer_oncology_flag(disease, explicit)
            for disease, explicit in zip(out["disease_term"], out["is_neoplasm"])
        ]
    else:
        out["is_neoplasm"] = out["disease_term"].map(infer_oncology_flag)
    for col in TABLE_SCHEMAS["death_mode_disease_pairs"]:
        if col not in out:
            out[col] = None
    out = out.drop_duplicates(subset=["death_mode", "disease_term"])
    return cast_for_schema(out, "death_mode_disease_pairs")


def normalize_disease_dictionary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "disease_term" not in df:
        return pd.DataFrame()
    out = df.copy()
    out["disease_term"] = out["disease_term"].map(normalize_disease)
    out = out[out["disease_term"] != ""].copy()
    if "is_neoplasm" in out:
        out["is_neoplasm"] = [
            infer_oncology_flag(disease, explicit)
            for disease, explicit in zip(out["disease_term"], out["is_neoplasm"])
        ]
    else:
        out["is_neoplasm"] = out["disease_term"].map(infer_oncology_flag)
    for col in ["disease_system", "tumor_family", "is_neoplasm"]:
        if col not in out:
            out[col] = None
    return out[["disease_term", "disease_system", "tumor_family", "is_neoplasm"]].drop_duplicates(
        subset=["disease_term"]
    )


def enrich_disease_attrs(df: pd.DataFrame, disease_dict: pd.DataFrame, oncology_col: str) -> pd.DataFrame:
    if df.empty or disease_dict.empty or "disease_term" not in df:
        return df
    out = df.merge(disease_dict, on="disease_term", how="left", suffixes=("", "_dict"))
    for target in ["disease_system", "tumor_family", oncology_col]:
        source = "is_neoplasm_dict" if target == oncology_col else f"{target}_dict"
        if target in out and source in out:
            current = out[target]
            out[target] = current.where(current.notna() & (current.astype(str) != ""), out[source])
    drop_cols = [col for col in out.columns if col.endswith("_dict")]
    if drop_cols:
        out = out.drop(columns=drop_cols)
    return out


def normalize_gene_mentions(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("gene_mentions")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    out["gene_symbol"] = out.get("gene_symbol", "").map(normalize_gene_symbol)
    death_mode = out.get("death_mode", pd.Series("", index=out.index))
    explicit = out.get("artifact_flag", pd.Series(None, index=out.index))
    out["artifact_flag"] = [
        as_bool(value) if as_bool(value) is not None else flag_gene_artifact(symbol, mode)
        for symbol, mode, value in zip(out["gene_symbol"], death_mode, explicit)
    ]
    out = out[(out["pmid"] != "") & (out["gene_symbol"] != "")].copy()
    for col in TABLE_SCHEMAS["gene_mentions"]:
        if col not in out:
            out[col] = None
    out = out.drop_duplicates(subset=["pmid", "gene_symbol", "source_field", "matched_text"])
    return cast_for_schema(out, "gene_mentions")


def normalize_gene_death_disease_articles(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("gene_death_disease_articles")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    out["gene_symbol"] = out.get("gene_symbol", "").map(normalize_gene_symbol)
    out["death_mode"] = out.get("death_mode", "").map(normalize_death_mode)
    out["disease_term"] = out.get("disease_term", "").map(normalize_disease)
    out = out[(out["pmid"] != "") & (out["gene_symbol"] != "")].copy()
    out = out[(out["death_mode"] == "") | out["death_mode"].isin(CANONICAL_DEATH_MODES)].copy()
    explicit = out.get("artifact_flag", pd.Series(None, index=out.index))
    out["artifact_flag"] = [
        as_bool(value) if as_bool(value) is not None else flag_gene_artifact(symbol, mode)
        for symbol, mode, value in zip(out["gene_symbol"], out["death_mode"], explicit)
    ]
    if "is_oncology" in out:
        out["is_oncology"] = [
            infer_oncology_flag(disease, explicit)
            for disease, explicit in zip(out["disease_term"], out["is_oncology"])
        ]
    else:
        out["is_oncology"] = out["disease_term"].map(infer_oncology_flag)
    for col in TABLE_SCHEMAS["gene_death_disease_articles"]:
        if col not in out:
            out[col] = None
    out = out.drop_duplicates(subset=["pmid", "gene_symbol", "death_mode", "disease_term", "context_sentence"])
    return cast_for_schema(out, "gene_death_disease_articles")


def normalize_drug_mentions(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("drug_mentions")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    if "drug_name" not in out:
        out["drug_name"] = out.get("normalized_drug_name", "")
    out["normalized_drug_name"] = out.get("normalized_drug_name", out["drug_name"]).map(normalize_drug_name)
    out = out[(out["pmid"] != "") & (out["normalized_drug_name"] != "")].copy()
    for col in TABLE_SCHEMAS["drug_mentions"]:
        if col not in out:
            out[col] = None
    out = _infer_drug_approval_flags(out)
    out = _prefer_approved_rows(out)
    out = out.drop_duplicates(subset=["pmid", "normalized_drug_name", "relation_sentence"])
    return cast_for_schema(out, "drug_mentions")


def normalize_drug_death_disease_articles(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return empty_frame("drug_death_disease_articles")
    out = df.copy()
    out["pmid"] = out.get("pmid", "").map(normalize_pmid)
    fallback_drug = out.get("drug_name", "")
    out["normalized_drug_name"] = out.get("normalized_drug_name", fallback_drug).map(normalize_drug_name)
    out["death_mode"] = out.get("death_mode", "").map(normalize_death_mode)
    out["disease_term"] = out.get("disease_term", "").map(normalize_disease)
    out = out[(out["pmid"] != "") & (out["normalized_drug_name"] != "")].copy()
    out = out[(out["death_mode"] == "") | out["death_mode"].isin(CANONICAL_DEATH_MODES)].copy()
    for col in TABLE_SCHEMAS["drug_death_disease_articles"]:
        if col not in out:
            out[col] = None
    out = _infer_drug_approval_flags(out)
    out = _prefer_approved_rows(out)
    out = out.drop_duplicates(
        subset=["pmid", "normalized_drug_name", "death_mode", "disease_term", "relation_sentence"]
    )
    return cast_for_schema(out, "drug_death_disease_articles")


def normalize_synonyms(df: pd.DataFrame, synonym_config: dict[str, Any]) -> pd.DataFrame:
    defaults: list[dict[str, str]] = []
    for canonical, synonyms in synonym_config.get("death_modes", {}).items():
        for synonym in synonyms:
            defaults.append(
                {"entity_type": "death_mode", "canonical_name": canonical, "synonym": str(synonym)}
            )
    out = pd.concat([pd.DataFrame(defaults), df], ignore_index=True, sort=False)
    if out.empty:
        return empty_frame("synonyms")
    for col in TABLE_SCHEMAS["synonyms"]:
        if col not in out:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str)
    out = out[(out["entity_type"] != "") & (out["canonical_name"] != "") & (out["synonym"] != "")]
    out = out.drop_duplicates(subset=["entity_type", "canonical_name", "synonym"])
    return cast_for_schema(out, "synonyms")


def derive_articles_from_pmids(tables: dict[str, pd.DataFrame], existing: pd.DataFrame) -> pd.DataFrame:
    pmids: set[str] = set(existing["pmid"].dropna().astype(str)) if not existing.empty else set()
    for df in tables.values():
        if "pmid" in df.columns:
            pmids.update(df["pmid"].dropna().astype(str))
    missing = sorted(pmid for pmid in pmids if pmid and pmid not in set(existing["pmid"].astype(str)))
    if not missing:
        return existing
    placeholders = pd.DataFrame({"pmid": missing})
    combined = pd.concat([existing, placeholders], ignore_index=True, sort=False)
    return normalize_articles(combined)


def derive_gene_mentions(gene_articles: pd.DataFrame, existing: pd.DataFrame) -> pd.DataFrame:
    if not existing.empty or gene_articles.empty:
        return existing
    cols = {
        "pmid": gene_articles["pmid"],
        "gene_symbol": gene_articles["gene_symbol"],
        "gene_alias": gene_articles["gene_symbol"],
        "source_field": gene_articles.get("source_field", ""),
        "matched_text": gene_articles["gene_symbol"],
        "context_sentence": gene_articles.get("context_sentence", ""),
        "artifact_flag": gene_articles.get("artifact_flag", False),
    }
    return normalize_gene_mentions(pd.DataFrame(cols))


def derive_drug_mentions(drug_articles: pd.DataFrame, existing: pd.DataFrame) -> pd.DataFrame:
    if not existing.empty or drug_articles.empty:
        return existing
    cols = {
        "pmid": drug_articles["pmid"],
        "drug_name": drug_articles["normalized_drug_name"],
        "normalized_drug_name": drug_articles["normalized_drug_name"],
        "reference_class_label": drug_articles.get("reference_class_label", ""),
        "approved_clinical_drug_flag": drug_articles.get("approved_clinical_drug_flag", False),
        "fda_approved_flag": drug_articles.get("fda_approved_flag", False),
        "source": "derived from drug_death_disease_articles",
        "relation_cue": drug_articles.get("relation_cue", ""),
        "relation_sentence": drug_articles.get("relation_sentence", ""),
        "negation_flag": drug_articles.get("negation_flag", False),
        "token_distance": drug_articles.get("token_distance", None),
        "relation_confidence": "",
    }
    return normalize_drug_mentions(pd.DataFrame(cols))


def enrich_articles_with_stage(articles: pd.DataFrame, article_stage: pd.DataFrame) -> pd.DataFrame:
    if articles.empty or article_stage.empty:
        return articles
    stage = article_stage.groupby("pmid", dropna=False).agg(
        inferred_oncology=("is_neoplasm", _first_nonempty),
        inferred_tumor_family=("tumor_family", _first_nonempty),
        inferred_evidence_stage=("evidence_stage", _first_nonempty),
    )
    out = articles.merge(stage.reset_index(), on="pmid", how="left")
    for target, source in [
        ("is_oncology", "inferred_oncology"),
        ("tumor_family", "inferred_tumor_family"),
        ("evidence_stage", "inferred_evidence_stage"),
    ]:
        out[target] = out[target].where(out[target].notna() & (out[target].astype(str) != ""), out[source])
    out = out.drop(columns=["inferred_oncology", "inferred_tumor_family", "inferred_evidence_stage"])
    return cast_for_schema(out, "articles")


def build_tables(
    project_dir: Path,
    raw_dir: Path,
    config_path: Path,
    synonym_path: Path,
) -> tuple[dict[str, pd.DataFrame], list[dict[str, Any]], list[str], dict[str, Any]]:
    config = load_yaml(config_path)
    synonym_config = load_yaml(synonym_path)
    use_demo = bool(config.get("database", {}).get("use_demo_if_no_raw", True))
    logical_configs = config.get("logical_tables", {})
    raw_tables: dict[str, pd.DataFrame] = {}
    source_records: list[dict[str, Any]] = []
    warnings: list[str] = []

    for logical_table, table_config in logical_configs.items():
        raw_df, records, load_warnings = load_logical_table(raw_dir, logical_table, table_config, use_demo)
        raw_tables[logical_table] = raw_df
        source_records.extend(records)
        warnings.extend(load_warnings)

    tables: dict[str, pd.DataFrame] = {}
    disease_dict = normalize_disease_dictionary(raw_tables.get("disease_dictionary", pd.DataFrame()))
    tables["death_mode_disease_articles"] = normalize_death_mode_disease_articles(
        raw_tables.get("death_mode_disease_articles", pd.DataFrame())
    )
    tables["death_mode_disease_articles"] = cast_for_schema(
        enrich_disease_attrs(tables["death_mode_disease_articles"], disease_dict, "is_neoplasm"),
        "death_mode_disease_articles",
    )
    tables["articles"] = normalize_articles(raw_tables.get("articles", pd.DataFrame()))
    tables["articles"] = derive_articles_from_pmids(tables, tables["articles"])
    tables["articles"] = enrich_articles_with_stage(tables["articles"], tables["death_mode_disease_articles"])
    tables["death_mode_disease_pairs"] = normalize_pairs(
        raw_tables.get("death_mode_disease_pairs", pd.DataFrame()),
        tables["death_mode_disease_articles"],
        tables["articles"],
    )
    tables["death_mode_disease_pairs"] = cast_for_schema(
        enrich_disease_attrs(tables["death_mode_disease_pairs"], disease_dict, "is_neoplasm"),
        "death_mode_disease_pairs",
    )
    tables["gene_death_disease_articles"] = normalize_gene_death_disease_articles(
        raw_tables.get("gene_death_disease_articles", pd.DataFrame())
    )
    tables["gene_death_disease_articles"] = cast_for_schema(
        enrich_disease_attrs(tables["gene_death_disease_articles"], disease_dict, "is_oncology"),
        "gene_death_disease_articles",
    )
    tables["gene_mentions"] = normalize_gene_mentions(raw_tables.get("gene_mentions", pd.DataFrame()))
    tables["gene_mentions"] = derive_gene_mentions(
        tables["gene_death_disease_articles"], tables["gene_mentions"]
    )
    tables["drug_death_disease_articles"] = normalize_drug_death_disease_articles(
        raw_tables.get("drug_death_disease_articles", pd.DataFrame())
    )
    tables["drug_mentions"] = normalize_drug_mentions(raw_tables.get("drug_mentions", pd.DataFrame()))
    tables["drug_mentions"] = derive_drug_mentions(
        tables["drug_death_disease_articles"], tables["drug_mentions"]
    )
    tables["articles"] = derive_articles_from_pmids(tables, tables["articles"])
    tables["articles"] = enrich_articles_with_stage(tables["articles"], tables["death_mode_disease_articles"])
    tables["synonyms"] = normalize_synonyms(raw_tables.get("synonyms", pd.DataFrame()), synonym_config)

    for table in table_names():
        if table not in tables:
            tables[table] = empty_frame(table)
        else:
            tables[table] = cast_for_schema(tables[table], table)

    return tables, source_records, warnings, config


def write_database(db_path: Path, tables: dict[str, pd.DataFrame], source_records: list[dict[str, Any]], validation_records: list[dict[str, str]]) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(db_path))
    try:
        create_tables(conn, drop=True)
        for table in table_names():
            insert_dataframe(conn, table, tables[table])
        insert_dataframe(conn, "build_sources", pd.DataFrame(source_records))
        insert_dataframe(conn, "validation_warnings", pd.DataFrame(validation_records))
        create_views(conn)
    finally:
        conn.close()


def export_static_tables(db_path: Path, exports_dir: Path, parquet_dir: Path) -> None:
    exports_dir.mkdir(parents=True, exist_ok=True)
    if parquet_dir.exists():
        shutil.rmtree(parquet_dir, ignore_errors=True)
    parquet_dir.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(db_path), read_only=True)
    try:
        for table in list(table_names()) + ["v_death_mode_summary", "v_disease_by_death_mode"]:
            df = conn.execute(f"SELECT * FROM {table}").fetchdf()
            df.to_csv(exports_dir / f"{table}.csv", index=False)
            parquet_path = str((parquet_dir / f"{table}.parquet").resolve()).replace("'", "''")
            conn.execute(f"COPY (SELECT * FROM {table}) TO '{parquet_path}' (FORMAT PARQUET)")
    finally:
        conn.close()


def build_database(
    project_dir: Path = PROJECT_ROOT,
    raw_dir: Path | None = None,
    config_path: Path | None = None,
    synonym_path: Path | None = None,
    db_path: Path | None = None,
    processed_dir: Path | None = None,
    exports_dir: Path | None = None,
    write_db: bool = True,
    export_outputs: bool = True,
) -> dict[str, Any]:
    project_dir = Path(project_dir)
    raw_dir = raw_dir or project_dir / "data" / "raw"
    config_path = config_path or project_dir / "config" / "schema_config.yaml"
    synonym_path = synonym_path or project_dir / "config" / "synonym_config.yaml"
    processed_dir = processed_dir or project_dir / "data" / "processed"
    exports_dir = exports_dir or project_dir / "data" / "exports"

    tables, source_records, warnings, config = build_tables(project_dir, raw_dir, config_path, synonym_path)
    year_cfg = config.get("year", {})
    validation_records, report = validate_tables(
        tables,
        source_records,
        warnings,
        year_min=int(year_cfg.get("min_expected", 2000)),
        year_max=int(year_cfg.get("max_expected", 2025)),
    )
    validation_path = processed_dir / "validation_report.md"
    write_validation_report(validation_path, report)

    if db_path is None:
        db_path = project_dir / config.get("database", {}).get("path", "data/processed/rcd_atlas.duckdb")
    if write_db:
        write_database(db_path, tables, source_records, validation_records)
        if export_outputs:
            parquet_dir = project_dir / config.get("database", {}).get(
                "parquet_dir", "data/processed/rcd_atlas.parquet"
            )
            export_static_tables(db_path, exports_dir, parquet_dir)

    return {
        "db_path": str(db_path),
        "validation_report": str(validation_path),
        "tables": {name: len(df) for name, df in tables.items()},
        "sources": source_records,
        "warnings": validation_records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the regulated cell-death atlas DuckDB database.")
    parser.add_argument("--project-dir", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--raw-dir", type=Path, default=None)
    parser.add_argument("--db-path", type=Path, default=None)
    parser.add_argument("--validate-only", action="store_true", help="Write validation report without DuckDB output.")
    parser.add_argument("--no-export", action="store_true", help="Skip CSV and Parquet export.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_database(
        project_dir=args.project_dir,
        raw_dir=args.raw_dir,
        db_path=args.db_path,
        write_db=not args.validate_only,
        export_outputs=not args.no_export,
    )
    print("Build complete")
    print(f"Database: {result['db_path']}")
    print(f"Validation report: {result['validation_report']}")
    print("Rows:")
    for table, rows in result["tables"].items():
        print(f"  {table}: {rows}")


if __name__ == "__main__":
    main()
