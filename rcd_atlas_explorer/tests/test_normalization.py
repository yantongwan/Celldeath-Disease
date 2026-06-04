from __future__ import annotations

from src.normalization import (
    flag_gene_artifact,
    make_pubmed_link,
    normalize_death_mode,
    normalize_gene_symbol,
    normalize_pmid,
)


def test_death_mode_normalization_maps_icd() -> None:
    assert normalize_death_mode("ICD") == "Immunogenic cell death"


def test_death_mode_normalization_maps_nets() -> None:
    assert normalize_death_mode("NETs") == "NETosis"


def test_pmids_are_strings() -> None:
    pmid = normalize_pmid("123456.0")
    assert pmid == "123456"
    assert isinstance(pmid, str)


def test_ascl4_artifact_flag() -> None:
    assert normalize_gene_symbol("ascl4") == "ASCL4"
    assert flag_gene_artifact("ASCL4", "Ferroptosis") is True


def test_pubmed_link() -> None:
    assert make_pubmed_link("1001") == "https://pubmed.ncbi.nlm.nih.gov/1001/"
