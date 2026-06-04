from __future__ import annotations

import tempfile
from pathlib import Path

import duckdb
import pandas as pd
import streamlit as st

from src import metrics
from src.normalization import CANONICAL_DEATH_MODES, normalize_death_mode, normalize_drug_name, normalize_gene_symbol
from src.queries import (
    get_articles_for_death_mode_disease,
    get_death_mode_disease_summary,
    global_article_search,
    query_drug_disease,
    query_drugs_for_death_mode_disease,
    query_gene_disease,
)
from src.schema import TABLE_SCHEMAS, create_tables, create_views
from src.ui_components import (
    DATA_DICTIONARY,
    apply_page_style,
    csv_bytes,
    download_table,
    load_ui_config,
    metric_row,
    pmid_text_area,
    safe_bar_chart,
    show_dataframe,
    summary_cards_from_articles,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "data" / "processed" / "rcd_atlas.duckdb"
PARQUET_DIR = PROJECT_ROOT / "data" / "processed" / "rcd_atlas.parquet"
RUNTIME_DB_PATH = Path(tempfile.gettempdir()) / "rcd_atlas_explorer.duckdb"


st.set_page_config(
    page_title="Interactive atlas of regulated cell death",
    layout="wide",
    initial_sidebar_state="collapsed",
)
apply_page_style()
UI_CONFIG = load_ui_config(PROJECT_ROOT)


@st.cache_data(show_spinner=False)
def read_sql(db_path: str, sql: str) -> pd.DataFrame:
    conn = duckdb.connect(db_path, read_only=True)
    try:
        return conn.execute(sql).fetchdf()
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def cached_overview(db_path: str) -> tuple[dict[str, int], pd.DataFrame, pd.DataFrame, tuple[int, int]]:
    conn = duckdb.connect(db_path, read_only=True)
    try:
        overview = metrics.overview_metrics(conn)
        summary = conn.execute("SELECT * FROM v_death_mode_summary ORDER BY total_unique_pmids DESC").fetchdf()
        counts = metrics.table_counts(conn)
        years = metrics.year_bounds(conn)
        return overview, summary, counts, years
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def distinct_values(db_path: str, table: str, column: str) -> list[str]:
    conn = duckdb.connect(db_path, read_only=True)
    try:
        return metrics.distinct_values(conn, table, column)
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def diseases_for_death_mode(db_path: str, death_mode: str) -> list[str]:
    conn = duckdb.connect(db_path, read_only=True)
    try:
        rows = conn.execute(
            """
            SELECT DISTINCT disease_term
            FROM death_mode_disease_articles
            WHERE death_mode = ?
              AND disease_term IS NOT NULL
              AND disease_term <> ''
            ORDER BY disease_term
            """,
            [death_mode],
        ).fetchall()
        return [str(row[0]) for row in rows]
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def diseases_for_query_context(
    db_path: str,
    source: str,
    death_modes: tuple[str, ...] = tuple(),
    disease_systems: tuple[str, ...] = tuple(),
    gene_symbol: str = "",
    gene_symbols: tuple[str, ...] = tuple(),
    drug_name: str = "",
    drug_names: tuple[str, ...] = tuple(),
    show_artifacts: bool = False,
) -> list[str]:
    modes = [normalize_death_mode(mode) for mode in death_modes if normalize_death_mode(mode)]
    systems = [system for system in disease_systems if system]
    where = ["d.disease_term IS NOT NULL", "d.disease_term <> ''"]
    params: list[str] = []

    if source == "gene":
        from_sql = "gene_death_disease_articles d"
        system_column = "d.disease_system"
        genes = []
        gene = normalize_gene_symbol(gene_symbol)
        if gene:
            genes.append(gene)
        genes.extend(normalize_gene_symbol(value) for value in gene_symbols)
        genes = sorted({value for value in genes if value})
        if genes:
            placeholders = ", ".join(["?"] * len(genes))
            where.append(f"d.gene_symbol IN ({placeholders})")
            params.extend(genes)
        if not show_artifacts:
            where.append("d.artifact_flag IS NOT TRUE")
    elif source == "drug":
        from_sql = """
            drug_death_disease_articles d
            LEFT JOIN death_mode_disease_pairs p
                ON p.death_mode = d.death_mode
                AND p.disease_term = d.disease_term
        """
        system_column = "p.disease_system"
        selected_drugs = [normalize_drug_name(value) for value in drug_names]
        selected_drugs = sorted({value for value in selected_drugs if value})
        if selected_drugs:
            placeholders = ", ".join(["?"] * len(selected_drugs))
            where.append(f"d.normalized_drug_name IN ({placeholders})")
            params.extend(selected_drugs)
        drug = normalize_drug_name(drug_name)
        if drug:
            where.append("LOWER(d.normalized_drug_name) LIKE ?")
            params.append(f"%{drug}%")
        where.append("d.negation_flag IS NOT TRUE")
    else:
        return []

    if modes:
        placeholders = ", ".join(["?"] * len(modes))
        where.append(f"d.death_mode IN ({placeholders})")
        params.extend(modes)
    if systems:
        placeholders = ", ".join(["?"] * len(systems))
        where.append(f"{system_column} IN ({placeholders})")
        params.extend(systems)

    conn = duckdb.connect(db_path, read_only=True)
    try:
        rows = conn.execute(
            f"""
            SELECT DISTINCT d.disease_term
            FROM {from_sql}
            WHERE {' AND '.join(where)}
            ORDER BY d.disease_term
            """,
            params,
        ).fetchall()
        return [str(row[0]) for row in rows]
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def cached_death_summary(
    db_path: str,
    death_mode: str,
    disease_search: str,
    min_pair_count: int,
    year_range: tuple[int, int],
    evidence_stages: tuple[str, ...],
    oncology_only: bool,
    sort_by: str,
) -> pd.DataFrame:
    return get_death_mode_disease_summary(
        death_mode=death_mode,
        disease_search=disease_search,
        min_pair_count=min_pair_count,
        year_range=year_range,
        evidence_stages=list(evidence_stages),
        oncology_only=oncology_only,
        sort_by=sort_by,
        db_path=db_path,
    )


@st.cache_data(show_spinner=False)
def cached_death_articles(
    db_path: str,
    death_mode: str,
    disease_term: str,
    year_range: tuple[int, int],
    evidence_stages: tuple[str, ...],
    oncology_only: bool,
) -> pd.DataFrame:
    return get_articles_for_death_mode_disease(
        death_mode,
        disease_term,
        {
            "year_range": year_range,
            "evidence_stages": list(evidence_stages),
            "oncology_only": oncology_only,
        },
        db_path=db_path,
    )


@st.cache_data(show_spinner=False)
def cached_gene_query(
    db_path: str,
    gene_symbols: tuple[str, ...],
    disease_search: str,
    death_modes: tuple[str, ...],
    evidence_stages: tuple[str, ...],
    source_fields: tuple[str, ...],
    disease_systems: tuple[str, ...],
    show_artifacts: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return query_gene_disease(
        gene_symbols=list(gene_symbols),
        disease_search=disease_search,
        death_modes=list(death_modes),
        filters={
            "evidence_stages": list(evidence_stages),
            "source_fields": list(source_fields),
            "disease_systems": list(disease_systems),
        },
        show_artifacts=show_artifacts,
        db_path=db_path,
    )


@st.cache_data(show_spinner=False)
def cached_drug_query(
    db_path: str,
    drug_names: tuple[str, ...],
    disease_search: str,
    death_modes: tuple[str, ...],
    approved_only: bool,
    reference_classes: tuple[str, ...],
    relation_cues: tuple[str, ...],
    disease_systems: tuple[str, ...],
    exclude_negated: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return query_drug_disease(
        drug_names=list(drug_names),
        disease_search=disease_search,
        death_modes=list(death_modes),
        approved_only=approved_only,
        reference_classes=list(reference_classes),
        relation_cues=list(relation_cues),
        disease_systems=list(disease_systems),
        exclude_negated=exclude_negated,
        db_path=db_path,
    )


@st.cache_data(show_spinner=False)
def cached_drugs_for_pair(
    db_path: str,
    death_mode: str,
    disease_term: str,
    approved_only: bool,
    reference_classes: tuple[str, ...],
    relation_cues: tuple[str, ...],
    disease_systems: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return query_drugs_for_death_mode_disease(
        death_mode=death_mode,
        disease_term=disease_term,
        approved_only=approved_only,
        reference_classes=list(reference_classes),
        relation_cues=list(relation_cues),
        disease_systems=list(disease_systems),
        db_path=db_path,
    )


@st.cache_data(show_spinner=False)
def cached_global_search(
    db_path: str,
    text: str,
    year_range: tuple[int, int],
    journals: tuple[str, ...],
    death_modes: tuple[str, ...],
    disease_systems: tuple[str, ...],
    evidence_stages: tuple[str, ...],
    oncology_choice: str,
    trial_supported: bool,
    gene_present: bool,
    drug_present: bool,
) -> pd.DataFrame:
    oncology = True if oncology_choice == "Oncology only" else False if oncology_choice == "Non-oncology only" else None
    return global_article_search(
        text=text,
        filters={
            "year_range": year_range,
            "journals": list(journals),
            "death_modes": list(death_modes),
            "disease_systems": list(disease_systems),
            "evidence_stages": list(evidence_stages),
            "oncology": oncology,
            "trial_supported": trial_supported,
            "gene_present": gene_present,
            "drug_present": drug_present,
            "limit": 5000,
        },
        db_path=db_path,
    )


@st.cache_resource(show_spinner="Preparing atlas database from bundled Parquet files...")
def materialize_runtime_database(parquet_dir: str, runtime_db_path: str) -> str:
    parquet_root = Path(parquet_dir)
    runtime_path = Path(runtime_db_path)
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    if runtime_path.exists():
        return str(runtime_path)

    conn = duckdb.connect(str(runtime_path))
    try:
        create_tables(conn, drop=True)
        for table, schema in TABLE_SCHEMAS.items():
            if table in {"build_sources", "validation_warnings"}:
                continue
            parquet_path = parquet_root / f"{table}.parquet"
            if not parquet_path.exists():
                continue
            columns = ", ".join(schema.keys())
            conn.execute(
                f"INSERT INTO {table} ({columns}) SELECT {columns} FROM read_parquet(?)",
                [str(parquet_path)],
            )
            rows_loaded = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            conn.execute(
                """
                INSERT INTO build_sources (logical_table, source_file, rows_loaded, column_mapping)
                VALUES (?, ?, ?, ?)
                """,
                [table, str(parquet_path.relative_to(PROJECT_ROOT)), int(rows_loaded), "bundled parquet"],
            )
        create_views(conn)
    finally:
        conn.close()
    return str(runtime_path)


def require_database() -> str:
    if DB_PATH.exists():
        return str(DB_PATH)
    if PARQUET_DIR.exists():
        return materialize_runtime_database(str(PARQUET_DIR), str(RUNTIME_DB_PATH))
    st.error("Database not found. Build it first with `python scripts/build_database.py`.")
    st.code("pip install -r requirements.txt\npython scripts/build_database.py\nstreamlit run app.py")
    return ""


st.title(UI_CONFIG.get("title", "Interactive atlas of regulated cell death across human disease"))
st.info(UI_CONFIG.get("warning", "This resource reports literature-level associations."))

DB = require_database()
if not DB:
    st.stop()

overview, death_summary, table_counts, year_bounds = cached_overview(DB)
death_modes = distinct_values(DB, "death_mode_disease_articles", "death_mode")
if not death_modes:
    death_modes = list(CANONICAL_DEATH_MODES)
disease_systems = distinct_values(DB, "death_mode_disease_articles", "disease_system")
gene_symbols = distinct_values(DB, "gene_death_disease_articles", "gene_symbol")
drug_names = distinct_values(DB, "drug_death_disease_articles", "normalized_drug_name")
evidence_stages = distinct_values(DB, "articles", "evidence_stage")
journals = distinct_values(DB, "articles", "journal")
source_fields = distinct_values(DB, "gene_mentions", "source_field")
reference_classes = distinct_values(DB, "drug_mentions", "reference_class_label")
relation_cues = distinct_values(DB, "drug_mentions", "relation_cue")

tabs = st.tabs(
    [
        "Atlas overview",
        "Death mode to diseases",
        "Gene + disease",
        "Drug / disease / death mode",
        "Article explorer",
        "Data dictionary",
    ]
)


with tabs[0]:
    metric_row(
        {
            "Unique PMIDs": overview.get("unique_pmids", 0),
            "Death modes": overview.get("death_modes", 0),
            "Diseases": overview.get("diseases", 0),
            "Pairs": overview.get("death_mode_disease_pairs", 0),
            "Genes": overview.get("genes", 0),
            "Drugs": overview.get("drugs", 0),
        }
    )
    left, right = st.columns(2)
    if not death_summary.empty:
        with left:
            st.subheader("PMIDs by death mode")
            safe_bar_chart(death_summary.set_index("death_mode")["total_unique_pmids"], "overview_pmid_chart")
        with right:
            st.subheader("Disease breadth by death mode")
            safe_bar_chart(death_summary.set_index("death_mode")["disease_breadth"], "overview_breadth_chart")
    st.subheader("Death-mode summary")
    show_dataframe(death_summary, "overview_summary")
    download_table(death_summary, "Download summary CSV", "death_mode_summary.csv", "overview_download")


with tabs[1]:
    c1, c2, c3, c4 = st.columns([1.1, 1.2, 1, 1])
    selected_mode = c1.selectbox("Death mode", death_modes, index=0, key="death_tab_mode")
    mode_diseases = diseases_for_death_mode(DB, selected_mode)
    selected_disease_filter = c2.selectbox(
        "Disease",
        ["All diseases"] + mode_diseases,
        key="death_tab_disease_filter",
    )
    disease_search = "" if selected_disease_filter == "All diseases" else selected_disease_filter
    selected_systems = c3.multiselect("Disease system", disease_systems, key="death_tab_disease_system")
    sort_by = c4.selectbox(
        "Sort by",
        ["pair_count", "disease_breadth", "latest_year", "literature_hub_score", "selectivity_score"],
        key="death_tab_sort_by",
    )
    c5, c6, c7, c8 = st.columns([1, 1.4, 1.2, 1])
    min_pair_count = c5.slider("Minimum pair count", 1, 50, 1, key="death_tab_min_pair_count")
    selected_years = c6.slider("Year range", year_bounds[0], year_bounds[1], year_bounds, key="death_tab_year_range")
    selected_evidence = c7.multiselect("Evidence stage", evidence_stages, key="death_tab_evidence")
    oncology_only = c8.toggle("Oncology only", value=False, key="death_tab_oncology_only")

    disease_df = cached_death_summary(
        DB,
        selected_mode,
        disease_search,
        min_pair_count,
        selected_years,
        tuple(selected_evidence),
        oncology_only,
        sort_by,
    )
    if selected_systems and "disease_system" in disease_df:
        disease_df = disease_df[disease_df["disease_system"].isin(selected_systems)]
    st.subheader("Disease summary")
    show_dataframe(disease_df, "death_disease_summary")
    d1, d2 = st.columns([1, 1])
    with d1:
        download_table(disease_df, "Download disease summary CSV", "death_mode_disease_summary.csv", "death_summary_download")
    with d2:
        pmid_text_area(pd.DataFrame(columns=["pmid"]), "death_summary_pmids")

    disease_options = disease_df["disease_term"].tolist() if "disease_term" in disease_df else []
    selected_disease = (
        st.selectbox("Select disease for article evidence", disease_options, key="death_tab_selected_disease")
        if disease_options
        else ""
    )
    if selected_disease:
        article_df = cached_death_articles(
            DB,
            selected_mode,
            selected_disease,
            selected_years,
            tuple(selected_evidence),
            oncology_only,
        )
        st.subheader("Article-level evidence")
        show_dataframe(article_df.head(500), "death_article_table")
        a1, a2 = st.columns([1, 1])
        with a1:
            download_table(article_df, "Download article evidence CSV", "death_mode_disease_articles.csv", "death_articles_download")
        with a2:
            pmid_text_area(article_df, "death_article_pmids")


with tabs[2]:
    c1, c2, c3, c4, c5 = st.columns([1.25, 1.15, 1.2, 1.35, 1])
    selected_gene_symbols = c1.multiselect("Gene symbol", gene_symbols, key="gene_tab_gene_symbols")
    gene_modes = c2.multiselect("Death mode", death_modes, key="gene_tab_death_modes")
    gene_systems = c3.multiselect("Disease system", disease_systems, key="gene_tab_disease_system")
    show_artifacts = c5.toggle("Show flagged artefacts", value=False, key="gene_tab_show_artifacts")
    gene_disease_options = diseases_for_query_context(
        DB,
        "gene",
        tuple(gene_modes),
        tuple(gene_systems),
        gene_symbols=tuple(selected_gene_symbols),
        show_artifacts=show_artifacts,
    )
    selected_gene_disease = c4.selectbox(
        "Disease",
        ["All diseases"] + gene_disease_options,
        key="gene_tab_disease_filter",
    )
    disease_input = "" if selected_gene_disease == "All diseases" else selected_gene_disease
    c6, c7 = st.columns(2)
    gene_evidence = c6.multiselect("Evidence stage", evidence_stages, key="gene_evidence")
    gene_sources = c7.multiselect("Source field", source_fields, key="gene_sources")
    gene_summary, gene_articles = cached_gene_query(
        DB,
        tuple(selected_gene_symbols),
        disease_input,
        tuple(gene_modes),
        tuple(gene_evidence),
        tuple(gene_sources),
        tuple(gene_systems),
        show_artifacts,
    )
    if (
        (not gene_summary.empty and "ASCL4" in set(gene_summary.get("gene_symbol", [])))
        or ("ASCL4" in {normalize_gene_symbol(value) for value in selected_gene_symbols})
    ):
        st.warning(
            "ASCL4 is flagged as a likely ACSL4 spelling/entity-normalization artefact in ferroptosis contexts and is excluded from biological summaries unless artefacts are enabled."
        )
    summary_cards_from_articles(gene_articles)
    if not gene_summary.empty:
        chart_left, chart_right = st.columns(2)
        with chart_left:
            st.subheader("Article counts by death mode")
            safe_bar_chart(gene_articles.groupby("death_mode")["pmid"].nunique(), "gene_death_mode_chart")
        with chart_right:
            st.subheader("Article counts by year")
            yearly = gene_articles.dropna(subset=["publication_year"]).groupby("publication_year")["pmid"].nunique()
            safe_bar_chart(yearly, "gene_year_chart")
    st.subheader("Gene-disease-death-mode summary")
    show_dataframe(gene_summary, "gene_summary_table")
    st.subheader("Article-level gene evidence")
    show_dataframe(gene_articles.head(500), "gene_articles_table")
    g1, g2 = st.columns([1, 1])
    with g1:
        download_table(gene_summary, "Download gene summary CSV", "gene_summary.csv", "gene_summary_download")
        download_table(gene_articles, "Download gene article CSV", "gene_articles.csv", "gene_articles_download")
    with g2:
        pmid_text_area(gene_articles, "gene_pmids")


with tabs[3]:
    st.warning(UI_CONFIG.get("candidate_drug_warning", "Drug rows are candidate literature associations."))
    mode = st.radio(
        "Query mode",
        ["Drug + disease", "Death mode + disease", "Drug only"],
        index=0,
        horizontal=True,
        key="drug_tab_query_mode",
    )
    f1, f2, f3, f4 = st.columns([1, 1.2, 1.2, 1.2])
    approved_only = f1.toggle("Approved clinical drug only", value=False, key="drug_tab_approved_only")
    selected_classes = f2.multiselect("Reference class", reference_classes, key="drug_tab_reference_class")
    selected_relation_cues = f3.multiselect("Relation cue", relation_cues, key="drug_tab_relation_cue")
    drug_systems = f4.multiselect("Disease system", disease_systems, key="drug_tab_disease_system")

    if mode == "Drug + disease":
        c1, c2, c3, c4 = st.columns([1.25, 1.2, 1.35, 1])
        selected_drug_names = c1.multiselect("Drug name", drug_names, key="drug_tab_a_drug_names")
        drug_modes = c2.multiselect("Death mode", death_modes, key="drug_mode_a")
        drug_disease_options = diseases_for_query_context(
            DB,
            "drug",
            tuple(drug_modes),
            tuple(drug_systems),
            drug_names=tuple(selected_drug_names),
        )
        selected_drug_disease = c3.selectbox(
            "Disease",
            ["All diseases"] + drug_disease_options,
            key="drug_tab_a_disease_filter",
        )
        drug_disease = "" if selected_drug_disease == "All diseases" else selected_drug_disease
        exclude_negated = c4.toggle("Exclude negated sentences", value=True, key="drug_tab_exclude_negated")
        drug_summary, drug_articles = cached_drug_query(
            DB,
            tuple(selected_drug_names),
            drug_disease,
            tuple(drug_modes),
            approved_only,
            tuple(selected_classes),
            tuple(selected_relation_cues),
            tuple(drug_systems),
            exclude_negated,
        )
    elif mode == "Death mode + disease":
        c1, c2 = st.columns([1, 1.4])
        pair_mode = c1.selectbox("Death mode", death_modes, key="drug_pair_mode")
        pair_disease_options = diseases_for_query_context(
            DB,
            "drug",
            (pair_mode,),
            tuple(drug_systems),
        )
        selected_pair_disease = c2.selectbox(
            "Disease",
            ["All diseases"] + pair_disease_options,
            key="drug_tab_b_disease_filter",
        )
        pair_disease = "" if selected_pair_disease == "All diseases" else selected_pair_disease
        drug_summary, drug_articles = cached_drugs_for_pair(
            DB,
            pair_mode,
            pair_disease,
            approved_only,
            tuple(selected_classes),
            tuple(selected_relation_cues),
            tuple(drug_systems),
        )
        drug_options = (
            sorted(set(drug_summary["normalized_drug_name"].dropna().astype(str)))
            if not drug_summary.empty and "normalized_drug_name" in drug_summary
            else []
        )
        selected_drug_for_pair = st.selectbox(
            "Select drug for article evidence",
            ["All drugs"] + drug_options,
            key="drug_tab_b_selected_drug",
        )
        if selected_drug_for_pair != "All drugs" and "normalized_drug_name" in drug_articles:
            drug_articles = drug_articles[drug_articles["normalized_drug_name"] == selected_drug_for_pair]
    else:
        selected_drug_names = st.multiselect("Drug name", drug_names, key="drug_tab_c_drug_names")
        drug_summary, drug_articles = cached_drug_query(
            DB,
            tuple(selected_drug_names),
            "",
            tuple(),
            approved_only,
            tuple(selected_classes),
            tuple(selected_relation_cues),
            tuple(drug_systems),
            True,
        )

    if not drug_articles.empty:
        col_a, col_b, col_c = st.columns(3)
        with col_a:
            st.subheader("Death-mode distribution")
            safe_bar_chart(drug_articles.groupby("death_mode")["pmid"].nunique(), "drug_death_mode_chart")
        with col_b:
            st.subheader("Disease distribution")
            safe_bar_chart(drug_articles.groupby("disease_term")["pmid"].nunique(), "drug_disease_chart")
        with col_c:
            st.subheader("Tumor-family distribution")
            safe_bar_chart(drug_articles.groupby("tumor_family")["pmid"].nunique(), "drug_tumor_family_chart")
    st.subheader("Drug summary")
    show_dataframe(drug_summary, "drug_summary_table")
    st.subheader("Article-level drug evidence")
    show_dataframe(drug_articles.head(500), "drug_articles_table")
    dr1, dr2 = st.columns([1, 1])
    with dr1:
        download_table(drug_summary, "Download drug summary CSV", "drug_summary.csv", "drug_summary_download")
        download_table(drug_articles, "Download drug article CSV", "drug_articles.csv", "drug_articles_download")
    with dr2:
        pmid_text_area(drug_articles, "drug_pmids")


with tabs[4]:
    c1, c2, c3 = st.columns([1.5, 1, 1])
    text = c1.text_input(
        "Search title, abstract, PMID, journal, death mode, disease, gene or drug",
        "",
        key="article_tab_text_search",
    )
    article_years = c2.slider("Year range", year_bounds[0], year_bounds[1], year_bounds, key="article_years")
    oncology_choice = c3.selectbox("Oncology", ["All", "Oncology only", "Non-oncology only"], key="article_tab_oncology")
    c4, c5, c6 = st.columns(3)
    selected_journals = c4.multiselect("Journal", journals, key="article_tab_journals")
    article_modes = c5.multiselect("Death mode", death_modes, key="article_modes")
    article_systems = c6.multiselect("Disease system", disease_systems, key="article_systems")
    c7, c8, c9 = st.columns(3)
    article_evidence = c7.multiselect("Evidence stage", evidence_stages, key="article_evidence")
    trial_supported = c8.toggle("Trial supported", value=False, key="article_tab_trial_supported")
    gene_present = c8.toggle("Gene present", value=False, key="article_tab_gene_present")
    drug_present = c9.toggle("Drug present", value=False, key="article_tab_drug_present")
    article_results = cached_global_search(
        DB,
        text,
        article_years,
        tuple(selected_journals),
        tuple(article_modes),
        tuple(article_systems),
        tuple(article_evidence),
        oncology_choice,
        trial_supported,
        gene_present,
        drug_present,
    )
    st.subheader("Articles")
    show_dataframe(article_results.head(500), "global_articles")
    with st.expander("Abstract snippets", expanded=False):
        for _, row in article_results.head(25).iterrows():
            st.markdown(f"**{row.get('pmid', '')} - {row.get('title', '')}**")
            st.write(row.get("abstract_snippet", ""))
    ar1, ar2 = st.columns([1, 1])
    with ar1:
        download_table(article_results, "Download selected rows CSV", "article_explorer_rows.csv", "article_download")
    with ar2:
        pmid_text_area(article_results, "article_pmids")


with tabs[5]:
    st.markdown(DATA_DICTIONARY)
    st.download_button(
        "Download data dictionary",
        data=DATA_DICTIONARY.encode("utf-8"),
        file_name="data_dictionary.md",
        mime="text/markdown",
    )
    st.subheader("Source files loaded")
    sources = read_sql(DB, "SELECT * FROM build_sources ORDER BY logical_table, source_file")
    show_dataframe(sources, "source_files_table", height=280)
    st.subheader("Rows per table")
    show_dataframe(table_counts, "table_counts", height=280)
    st.subheader("Validation warnings")
    warnings = read_sql(DB, "SELECT * FROM validation_warnings ORDER BY severity, message")
    show_dataframe(warnings, "validation_warnings", height=320)
    validation_report = PROJECT_ROOT / "data" / "processed" / "validation_report.md"
    if validation_report.exists():
        st.download_button(
            "Download validation report",
            data=validation_report.read_bytes(),
            file_name="validation_report.md",
            mime="text/markdown",
        )
