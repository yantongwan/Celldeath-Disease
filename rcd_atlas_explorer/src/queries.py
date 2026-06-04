from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import duckdb
import pandas as pd

from .normalization import normalize_death_mode, normalize_drug_name, normalize_gene_symbol


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = PROJECT_ROOT / "data" / "processed" / "rcd_atlas.duckdb"


DEATH_SUMMARY_COLUMNS = [
    "disease_term",
    "disease_system",
    "pair_count",
    "unique_pmids",
    "first_year",
    "latest_year",
    "literature_hub_score",
    "selectivity_score",
    "number_of_articles_with_genes",
    "number_of_articles_with_drugs",
    "trial_supported_count",
]

ARTICLE_COLUMNS = [
    "pmid",
    "title",
    "journal",
    "publication_year",
    "death_mode",
    "disease_term",
    "evidence_stage",
    "publication_types",
    "country",
    "pubmed_link",
    "doi_link",
    "abstract_snippet",
]

GENE_SUMMARY_COLUMNS = [
    "gene_symbol",
    "disease_term",
    "disease_system",
    "death_mode",
    "unique_pmids",
    "first_year",
    "latest_year",
    "evidence_stage_distribution",
    "oncology_fraction",
    "artifact_flag",
]

GENE_ARTICLE_COLUMNS = [
    "pmid",
    "title",
    "journal",
    "publication_year",
    "gene_symbol",
    "death_mode",
    "disease_term",
    "disease_system",
    "source_field",
    "context_sentence",
    "evidence_stage",
    "artifact_flag",
    "pubmed_link",
]

DRUG_SUMMARY_COLUMNS = [
    "normalized_drug_name",
    "disease_term",
    "disease_system",
    "death_mode",
    "reference_class_label",
    "approved_clinical_drug_flag",
    "unique_pmids",
    "first_year",
    "latest_year",
    "relation_cues",
    "top_relation_sentence",
]

DRUG_ARTICLE_COLUMNS = [
    "pmid",
    "title",
    "journal",
    "publication_year",
    "normalized_drug_name",
    "death_mode",
    "disease_term",
    "disease_system",
    "tumor_family",
    "reference_class_label",
    "approved_clinical_drug_flag",
    "relation_cue",
    "relation_sentence",
    "negation_flag",
    "evidence_stage",
    "pubmed_link",
]

GLOBAL_ARTICLE_COLUMNS = [
    "pmid",
    "title",
    "journal",
    "publication_year",
    "publication_types",
    "country",
    "evidence_stage",
    "is_oncology",
    "publicationtype_supported_trial",
    "death_modes",
    "diseases",
    "disease_systems",
    "genes",
    "drugs",
    "pubmed_link",
    "doi_link",
    "abstract_snippet",
]


def _empty(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def _connect(
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> tuple[duckdb.DuckDBPyConnection, bool]:
    if conn is not None:
        return conn, False
    path = Path(db_path) if db_path else DEFAULT_DB_PATH
    return duckdb.connect(str(path), read_only=True), True


def _fetch_df(
    sql: str,
    params: list[Any] | None = None,
    columns: Iterable[str] | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> pd.DataFrame:
    handle, should_close = _connect(conn=conn, db_path=db_path)
    try:
        return handle.execute(sql, params or []).fetchdf()
    except Exception:
        if columns is not None:
            return _empty(columns)
        raise
    finally:
        if should_close:
            handle.close()


def _append_in(where: list[str], params: list[Any], column: str, values: list[str] | None) -> None:
    cleaned = [value for value in (values or []) if value]
    if not cleaned:
        return
    placeholders = ", ".join(["?"] * len(cleaned))
    where.append(f"{column} IN ({placeholders})")
    params.extend(cleaned)


def _like(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip().lower()
    return f"%{value}%" if value else None


def _pubmed_expr(alias: str = "a") -> str:
    return f"CASE WHEN {alias}.pmid IS NOT NULL AND {alias}.pmid <> '' THEN 'https://pubmed.ncbi.nlm.nih.gov/' || {alias}.pmid || '/' ELSE '' END"


def _doi_expr(alias: str = "a") -> str:
    return (
        f"CASE WHEN {alias}.doi IS NULL OR {alias}.doi = '' THEN '' "
        f"WHEN STARTS_WITH({alias}.doi, 'http://') OR STARTS_WITH({alias}.doi, 'https://') THEN {alias}.doi "
        f"ELSE 'https://doi.org/' || {alias}.doi END"
    )


def _snippet_expr(alias: str = "a") -> str:
    return f"SUBSTR(COALESCE({alias}.abstract, ''), 1, 280)"


def get_death_mode_disease_summary(
    death_mode: str,
    disease_search: str | None = None,
    min_pair_count: int = 1,
    year_range: tuple[int, int] | None = None,
    evidence_stages: list[str] | None = None,
    oncology_only: bool = False,
    sort_by: str = "pair_count",
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> pd.DataFrame:
    mode = normalize_death_mode(death_mode)
    if not mode:
        return _empty(DEATH_SUMMARY_COLUMNS)
    where = ["d.death_mode = ?"]
    params: list[Any] = [mode]
    search = _like(disease_search)
    if search:
        where.append("LOWER(d.disease_term) LIKE ?")
        params.append(search)
    if year_range:
        where.append("a.publication_year BETWEEN ? AND ?")
        params.extend([int(year_range[0]), int(year_range[1])])
    _append_in(where, params, "COALESCE(d.evidence_stage, a.evidence_stage)", evidence_stages)
    if oncology_only:
        where.append("COALESCE(a.is_oncology, d.is_neoplasm) IS TRUE")
    sort_map = {
        "pair_count": "pair_count",
        "disease_breadth": "unique_pmids",
        "latest_year": "latest_year",
        "literature_hub_score": "literature_hub_score",
        "selectivity_score": "selectivity_score",
    }
    order_col = sort_map.get(sort_by, "pair_count")
    sql = f"""
        WITH base AS (
            SELECT
                d.pmid,
                d.death_mode,
                d.disease_term,
                d.disease_system,
                d.is_neoplasm,
                COALESCE(d.evidence_stage, a.evidence_stage) AS evidence_stage,
                a.publication_year,
                a.publicationtype_supported_trial
            FROM death_mode_disease_articles d
            LEFT JOIN articles a ON a.pmid = d.pmid
            WHERE {' AND '.join(where)}
        )
        SELECT
            b.disease_term,
            COALESCE(p.disease_system, ANY_VALUE(b.disease_system), '') AS disease_system,
            COALESCE(p.pair_count, COUNT(DISTINCT b.pmid)) AS pair_count,
            COALESCE(p.unique_pmids, COUNT(DISTINCT b.pmid)) AS unique_pmids,
            COALESCE(p.first_year, MIN(b.publication_year)) AS first_year,
            COALESCE(p.latest_year, MAX(b.publication_year)) AS latest_year,
            p.literature_hub_score,
            p.selectivity_score,
            COUNT(DISTINCT CASE WHEN g.pmid IS NOT NULL THEN b.pmid ELSE NULL END)
                AS number_of_articles_with_genes,
            COUNT(DISTINCT CASE WHEN dr.pmid IS NOT NULL THEN b.pmid ELSE NULL END)
                AS number_of_articles_with_drugs,
            COUNT(DISTINCT CASE WHEN b.publicationtype_supported_trial IS TRUE THEN b.pmid ELSE NULL END)
                AS trial_supported_count
        FROM base b
        LEFT JOIN death_mode_disease_pairs p
            ON p.death_mode = b.death_mode AND p.disease_term = b.disease_term
        LEFT JOIN gene_death_disease_articles g
            ON g.pmid = b.pmid
            AND g.death_mode = b.death_mode
            AND g.disease_term = b.disease_term
            AND g.artifact_flag IS NOT TRUE
        LEFT JOIN drug_death_disease_articles dr
            ON dr.pmid = b.pmid
            AND dr.death_mode = b.death_mode
            AND dr.disease_term = b.disease_term
        GROUP BY
            b.disease_term,
            p.disease_system,
            p.pair_count,
            p.unique_pmids,
            p.first_year,
            p.latest_year,
            p.literature_hub_score,
            p.selectivity_score
        HAVING COALESCE(p.pair_count, COUNT(DISTINCT b.pmid)) >= ?
        ORDER BY {order_col} DESC NULLS LAST, b.disease_term
    """
    params.append(int(min_pair_count))
    return _fetch_df(sql, params, DEATH_SUMMARY_COLUMNS, conn=conn, db_path=db_path)


def get_articles_for_death_mode_disease(
    death_mode: str,
    disease_term: str,
    filters: dict[str, Any] | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> pd.DataFrame:
    filters = filters or {}
    mode = normalize_death_mode(death_mode)
    disease = disease_term.strip() if disease_term else ""
    if not mode or not disease:
        return _empty(ARTICLE_COLUMNS)
    where = ["d.death_mode = ?", "d.disease_term = ?"]
    params: list[Any] = [mode, disease]
    year_range = filters.get("year_range")
    if year_range:
        where.append("a.publication_year BETWEEN ? AND ?")
        params.extend([int(year_range[0]), int(year_range[1])])
    _append_in(where, params, "COALESCE(d.evidence_stage, a.evidence_stage)", filters.get("evidence_stages"))
    if filters.get("oncology_only"):
        where.append("COALESCE(a.is_oncology, d.is_neoplasm) IS TRUE")
    sql = f"""
        SELECT DISTINCT
            d.pmid,
            COALESCE(a.title, '') AS title,
            COALESCE(a.journal, '') AS journal,
            a.publication_year,
            d.death_mode,
            d.disease_term,
            COALESCE(d.evidence_stage, a.evidence_stage, '') AS evidence_stage,
            COALESCE(a.publication_types, '') AS publication_types,
            COALESCE(a.country, '') AS country,
            {_pubmed_expr('a')} AS pubmed_link,
            {_doi_expr('a')} AS doi_link,
            {_snippet_expr('a')} AS abstract_snippet
        FROM death_mode_disease_articles d
        LEFT JOIN articles a ON a.pmid = d.pmid
        WHERE {' AND '.join(where)}
        ORDER BY a.publication_year DESC NULLS LAST, d.pmid
    """
    return _fetch_df(sql, params, ARTICLE_COLUMNS, conn=conn, db_path=db_path)


def query_gene_disease(
    gene_symbol: str | None = None,
    disease_search: str | None = None,
    death_modes: list[str] | None = None,
    filters: dict[str, Any] | None = None,
    gene_symbols: list[str] | None = None,
    show_artifacts: bool = False,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    filters = filters or {}
    where: list[str] = []
    params: list[Any] = []
    has_user_filter = False
    genes = []
    gene = normalize_gene_symbol(gene_symbol)
    if gene:
        genes.append(gene)
    genes.extend(normalize_gene_symbol(value) for value in (gene_symbols or []))
    genes = sorted({value for value in genes if value})
    _append_in(where, params, "g.gene_symbol", genes)
    if genes:
        has_user_filter = True
    search = _like(disease_search)
    if search:
        where.append("LOWER(g.disease_term) LIKE ?")
        params.append(search)
        has_user_filter = True
    modes = [normalize_death_mode(mode) for mode in (death_modes or []) if normalize_death_mode(mode)]
    _append_in(where, params, "g.death_mode", modes)
    has_user_filter = has_user_filter or bool(modes)
    evidence_stages = filters.get("evidence_stages")
    source_fields = filters.get("source_fields")
    disease_systems = [value for value in (filters.get("disease_systems") or []) if value]
    _append_in(where, params, "g.evidence_stage", evidence_stages)
    _append_in(where, params, "g.source_field", source_fields)
    _append_in(where, params, "g.disease_system", disease_systems)
    has_user_filter = has_user_filter or bool(evidence_stages) or bool(source_fields) or bool(disease_systems)
    if not has_user_filter:
        return _empty(GENE_SUMMARY_COLUMNS), _empty(GENE_ARTICLE_COLUMNS)
    if not show_artifacts:
        where.append("g.artifact_flag IS NOT TRUE")
    where_sql = " AND ".join(where)
    base_sql = f"""
        FROM gene_death_disease_articles g
        LEFT JOIN articles a ON a.pmid = g.pmid
        WHERE {where_sql}
    """
    summary_sql = f"""
        WITH base AS (
            SELECT
                g.pmid,
                g.gene_symbol,
                g.disease_term,
                COALESCE(g.disease_system, '') AS disease_system,
                g.death_mode,
                COALESCE(g.evidence_stage, a.evidence_stage, 'missing') AS evidence_stage,
                g.is_oncology,
                g.artifact_flag,
                a.publication_year
            {base_sql}
        ),
        stage_counts AS (
            SELECT
                gene_symbol,
                disease_term,
                disease_system,
                death_mode,
                evidence_stage,
                COUNT(DISTINCT pmid) AS n
            FROM base
            GROUP BY gene_symbol, disease_term, disease_system, death_mode, evidence_stage
        ),
        stage_text AS (
            SELECT
                gene_symbol,
                disease_term,
                disease_system,
                death_mode,
                STRING_AGG(evidence_stage || '=' || CAST(n AS VARCHAR), '; ' ORDER BY evidence_stage)
                    AS evidence_stage_distribution
            FROM stage_counts
            GROUP BY gene_symbol, disease_term, disease_system, death_mode
        )
        SELECT
            b.gene_symbol,
            b.disease_term,
            b.disease_system,
            b.death_mode,
            COUNT(DISTINCT b.pmid) AS unique_pmids,
            MIN(b.publication_year) AS first_year,
            MAX(b.publication_year) AS latest_year,
            COALESCE(s.evidence_stage_distribution, '') AS evidence_stage_distribution,
            AVG(CASE WHEN b.is_oncology IS TRUE THEN 1.0 ELSE 0.0 END) AS oncology_fraction,
            BOOL_OR(b.artifact_flag) AS artifact_flag
        FROM base b
        LEFT JOIN stage_text s
            ON s.gene_symbol = b.gene_symbol
            AND s.disease_term = b.disease_term
            AND s.disease_system = b.disease_system
            AND s.death_mode = b.death_mode
        GROUP BY b.gene_symbol, b.disease_term, b.disease_system, b.death_mode, s.evidence_stage_distribution
        ORDER BY unique_pmids DESC, b.gene_symbol, b.disease_term
    """
    article_sql = f"""
        SELECT DISTINCT
            g.pmid,
            COALESCE(a.title, '') AS title,
            COALESCE(a.journal, '') AS journal,
            a.publication_year,
            g.gene_symbol,
            g.death_mode,
            g.disease_term,
            COALESCE(g.disease_system, '') AS disease_system,
            COALESCE(g.source_field, '') AS source_field,
            COALESCE(g.context_sentence, '') AS context_sentence,
            COALESCE(g.evidence_stage, a.evidence_stage, '') AS evidence_stage,
            g.artifact_flag,
            {_pubmed_expr('a')} AS pubmed_link
        {base_sql}
        ORDER BY a.publication_year DESC NULLS LAST, g.pmid
    """
    summary = _fetch_df(summary_sql, params, GENE_SUMMARY_COLUMNS, conn=conn, db_path=db_path)
    articles = _fetch_df(article_sql, params, GENE_ARTICLE_COLUMNS, conn=conn, db_path=db_path)
    return summary, articles


def query_drug_disease(
    drug_name: str | None = None,
    disease_search: str | None = None,
    death_modes: list[str] | None = None,
    approved_only: bool = False,
    reference_classes: list[str] | None = None,
    relation_cues: list[str] | None = None,
    exclude_negated: bool = True,
    disease_systems: list[str] | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
    drug_names: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    where: list[str] = []
    params: list[Any] = []
    has_user_filter = False
    selected_drugs = [normalize_drug_name(value) for value in (drug_names or [])]
    selected_drugs = sorted({value for value in selected_drugs if value})
    _append_in(where, params, "d.normalized_drug_name", selected_drugs)
    if selected_drugs:
        has_user_filter = True
    drug = normalize_drug_name(drug_name)
    if drug:
        where.append("LOWER(d.normalized_drug_name) LIKE ?")
        params.append(f"%{drug}%")
        has_user_filter = True
    search = _like(disease_search)
    if search:
        where.append("LOWER(d.disease_term) LIKE ?")
        params.append(search)
        has_user_filter = True
    modes = [normalize_death_mode(mode) for mode in (death_modes or []) if normalize_death_mode(mode)]
    _append_in(where, params, "d.death_mode", modes)
    has_user_filter = has_user_filter or bool(modes)
    if approved_only:
        where.append("(d.approved_clinical_drug_flag IS TRUE OR d.fda_approved_flag IS TRUE)")
        has_user_filter = True
    _append_in(where, params, "d.reference_class_label", reference_classes)
    _append_in(where, params, "d.relation_cue", relation_cues)
    systems = [value for value in (disease_systems or []) if value]
    _append_in(where, params, "p.disease_system", systems)
    has_user_filter = has_user_filter or bool(reference_classes) or bool(relation_cues) or bool(systems)
    if not has_user_filter:
        return _empty(DRUG_SUMMARY_COLUMNS), _empty(DRUG_ARTICLE_COLUMNS)
    if exclude_negated:
        where.append("d.negation_flag IS NOT TRUE")
    where_sql = " AND ".join(where)
    base_sql = f"""
        FROM drug_death_disease_articles d
        LEFT JOIN articles a ON a.pmid = d.pmid
        LEFT JOIN death_mode_disease_pairs p
            ON p.death_mode = d.death_mode
            AND p.disease_term = d.disease_term
        WHERE {where_sql}
    """
    summary_sql = f"""
        SELECT
            d.normalized_drug_name,
            d.disease_term,
            COALESCE(p.disease_system, '') AS disease_system,
            d.death_mode,
            STRING_AGG(DISTINCT COALESCE(d.reference_class_label, ''), '; ') AS reference_class_label,
            BOOL_OR(d.approved_clinical_drug_flag) AS approved_clinical_drug_flag,
            COUNT(DISTINCT d.pmid) AS unique_pmids,
            MIN(a.publication_year) AS first_year,
            MAX(a.publication_year) AS latest_year,
            STRING_AGG(DISTINCT COALESCE(d.relation_cue, ''), '; ') AS relation_cues,
            ANY_VALUE(COALESCE(d.relation_sentence, '')) AS top_relation_sentence
        {base_sql}
        GROUP BY d.normalized_drug_name, d.disease_term, COALESCE(p.disease_system, ''), d.death_mode
        ORDER BY unique_pmids DESC, d.normalized_drug_name, d.disease_term
    """
    article_sql = f"""
        SELECT DISTINCT
            d.pmid,
            COALESCE(a.title, '') AS title,
            COALESCE(a.journal, '') AS journal,
            a.publication_year,
            d.normalized_drug_name,
            d.death_mode,
            d.disease_term,
            COALESCE(p.disease_system, '') AS disease_system,
            COALESCE(d.tumor_family, '') AS tumor_family,
            COALESCE(d.reference_class_label, '') AS reference_class_label,
            d.approved_clinical_drug_flag,
            COALESCE(d.relation_cue, '') AS relation_cue,
            COALESCE(d.relation_sentence, '') AS relation_sentence,
            d.negation_flag,
            COALESCE(d.evidence_stage, a.evidence_stage, '') AS evidence_stage,
            {_pubmed_expr('a')} AS pubmed_link
        {base_sql}
        ORDER BY a.publication_year DESC NULLS LAST, d.pmid
    """
    summary = _fetch_df(summary_sql, params, DRUG_SUMMARY_COLUMNS, conn=conn, db_path=db_path)
    articles = _fetch_df(article_sql, params, DRUG_ARTICLE_COLUMNS, conn=conn, db_path=db_path)
    return summary, articles


def query_drugs_for_death_mode_disease(
    death_mode: str,
    disease_term: str,
    approved_only: bool = False,
    reference_classes: list[str] | None = None,
    relation_cues: list[str] | None = None,
    disease_systems: list[str] | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return query_drug_disease(
        drug_name=None,
        disease_search=disease_term,
        death_modes=[death_mode],
        approved_only=approved_only,
        reference_classes=reference_classes,
        relation_cues=relation_cues,
        exclude_negated=True,
        disease_systems=disease_systems,
        conn=conn,
        db_path=db_path,
    )


def global_article_search(
    text: str | None = None,
    filters: dict[str, Any] | None = None,
    conn: duckdb.DuckDBPyConnection | None = None,
    db_path: str | Path | None = None,
) -> pd.DataFrame:
    filters = filters or {}
    where: list[str] = []
    params: list[Any] = []
    search = _like(text)
    if search:
        where.append(
            "LOWER(COALESCE(a.title, '') || ' ' || COALESCE(a.abstract, '') || ' ' || "
            "a.pmid || ' ' || COALESCE(a.journal, '') || ' ' || COALESCE(death.death_modes, '') || ' ' || "
            "COALESCE(death.diseases, '') || ' ' || COALESCE(gene.genes, '') || ' ' || "
            "COALESCE(drug.drugs, '')) LIKE ?"
        )
        params.append(search)
    year_range = filters.get("year_range")
    if year_range:
        where.append("a.publication_year BETWEEN ? AND ?")
        params.extend([int(year_range[0]), int(year_range[1])])
    _append_in(where, params, "a.journal", filters.get("journals"))
    _append_in(where, params, "a.evidence_stage", filters.get("evidence_stages"))
    oncology = filters.get("oncology")
    if oncology is True:
        where.append("a.is_oncology IS TRUE")
    elif oncology is False:
        where.append("a.is_oncology IS NOT TRUE")
    if filters.get("trial_supported"):
        where.append("a.publicationtype_supported_trial IS TRUE")
    modes = [normalize_death_mode(mode) for mode in filters.get("death_modes", []) if normalize_death_mode(mode)]
    if modes:
        placeholders = ", ".join(["?"] * len(modes))
        where.append(
            f"EXISTS (SELECT 1 FROM death_mode_disease_articles f WHERE f.pmid = a.pmid AND f.death_mode IN ({placeholders}))"
        )
        params.extend(modes)
    systems = filters.get("disease_systems") or []
    if systems:
        placeholders = ", ".join(["?"] * len(systems))
        where.append(
            f"EXISTS (SELECT 1 FROM death_mode_disease_articles f WHERE f.pmid = a.pmid AND f.disease_system IN ({placeholders}))"
        )
        params.extend(systems)
    if filters.get("gene_present"):
        where.append("gene.genes IS NOT NULL AND gene.genes <> ''")
    if filters.get("drug_present"):
        where.append("drug.drugs IS NOT NULL AND drug.drugs <> ''")
    where_sql = f"WHERE {' AND '.join(where)}" if where else ""
    limit = int(filters.get("limit", 5000))
    sql = f"""
        WITH death AS (
            SELECT
                pmid,
                STRING_AGG(DISTINCT death_mode, '; ' ORDER BY death_mode) AS death_modes,
                STRING_AGG(DISTINCT disease_term, '; ' ORDER BY disease_term) AS diseases,
                STRING_AGG(DISTINCT COALESCE(disease_system, ''), '; ' ORDER BY COALESCE(disease_system, ''))
                    AS disease_systems
            FROM death_mode_disease_articles
            GROUP BY pmid
        ),
        gene AS (
            SELECT
                pmid,
                STRING_AGG(DISTINCT gene_symbol, '; ' ORDER BY gene_symbol) AS genes
            FROM gene_death_disease_articles
            WHERE artifact_flag IS NOT TRUE
            GROUP BY pmid
        ),
        drug AS (
            SELECT
                pmid,
                STRING_AGG(DISTINCT normalized_drug_name, '; ' ORDER BY normalized_drug_name) AS drugs
            FROM drug_death_disease_articles
            WHERE negation_flag IS NOT TRUE
            GROUP BY pmid
        )
        SELECT
            a.pmid,
            COALESCE(a.title, '') AS title,
            COALESCE(a.journal, '') AS journal,
            a.publication_year,
            COALESCE(a.publication_types, '') AS publication_types,
            COALESCE(a.country, '') AS country,
            COALESCE(a.evidence_stage, '') AS evidence_stage,
            a.is_oncology,
            a.publicationtype_supported_trial,
            COALESCE(death.death_modes, '') AS death_modes,
            COALESCE(death.diseases, '') AS diseases,
            COALESCE(death.disease_systems, '') AS disease_systems,
            COALESCE(gene.genes, '') AS genes,
            COALESCE(drug.drugs, '') AS drugs,
            {_pubmed_expr('a')} AS pubmed_link,
            {_doi_expr('a')} AS doi_link,
            {_snippet_expr('a')} AS abstract_snippet
        FROM articles a
        LEFT JOIN death ON death.pmid = a.pmid
        LEFT JOIN gene ON gene.pmid = a.pmid
        LEFT JOIN drug ON drug.pmid = a.pmid
        {where_sql}
        ORDER BY a.publication_year DESC NULLS LAST, a.pmid
        LIMIT {limit}
    """
    return _fetch_df(sql, params, GLOBAL_ARTICLE_COLUMNS, conn=conn, db_path=db_path)
