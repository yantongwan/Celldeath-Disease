from __future__ import annotations

import fnmatch
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


SUPPORTED_SUFFIXES = {".csv", ".tsv", ".txt", ".xlsx", ".xls"}


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def normalize_header(value: str) -> str:
    value = value.strip().replace("\ufeff", "")
    value = re.sub(r"(?<!^)(?=[A-Z][a-z])", "_", value)
    value = re.sub(r"[^A-Za-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_").lower()


def is_hidden_or_sidecar(path: Path) -> bool:
    return any(part.startswith(".") or part.startswith("._") for part in path.parts)


def iter_raw_files(raw_dir: Path, include_demo: bool) -> list[Path]:
    if not raw_dir.exists():
        return []
    files: list[Path] = []
    for path in raw_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in SUPPORTED_SUFFIXES:
            continue
        rel_parts = path.relative_to(raw_dir).parts
        if is_hidden_or_sidecar(path):
            continue
        if not include_demo and "demo" in rel_parts:
            continue
        files.append(path)
    return sorted(files)


def choose_raw_pool(raw_dir: Path, use_demo_if_no_raw: bool = True) -> list[Path]:
    non_demo = iter_raw_files(raw_dir, include_demo=False)
    if non_demo:
        return non_demo
    if use_demo_if_no_raw:
        return iter_raw_files(raw_dir, include_demo=True)
    return []


def read_table_file(path: Path, nrows: int | None = None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False, nrows=nrows)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, nrows=nrows)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str, nrows=nrows).fillna("")
    raise ValueError(f"Unsupported input file type: {path}")


def build_alias_lookup(table_config: dict[str, Any]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for canonical, aliases in table_config.get("columns", {}).items():
        lookup[normalize_header(canonical)] = canonical
        for alias in aliases or []:
            lookup[normalize_header(str(alias))] = canonical
    return lookup


def standardize_columns(
    df: pd.DataFrame, table_config: dict[str, Any]
) -> tuple[pd.DataFrame, dict[str, str]]:
    alias_lookup = build_alias_lookup(table_config)
    exact_lookup = {
        normalize_header(canonical): canonical for canonical in table_config.get("columns", {})
    }
    out = pd.DataFrame(index=df.index)
    mapping: dict[str, str] = {}

    for source_col in df.columns:
        canonical = exact_lookup.get(normalize_header(str(source_col)))
        if canonical and canonical not in out.columns:
            out[canonical] = df[source_col]
            mapping[canonical] = str(source_col)

    for source_col in df.columns:
        canonical = alias_lookup.get(normalize_header(str(source_col)))
        if canonical and canonical not in out.columns:
            out[canonical] = df[source_col]
            mapping[canonical] = str(source_col)
    return out, mapping


def has_required_columns(mapping: dict[str, str], table_config: dict[str, Any]) -> bool:
    required = table_config.get("required_columns", [])
    return all(col in mapping for col in required)


def detect_source_files(
    raw_dir: Path, logical_table: str, table_config: dict[str, Any], use_demo_if_no_raw: bool
) -> tuple[list[Path], list[str]]:
    warnings: list[str] = []
    pool = choose_raw_pool(raw_dir, use_demo_if_no_raw=use_demo_if_no_raw)
    patterns = [pattern.lower() for pattern in table_config.get("filename_patterns", [])]
    matched: list[Path] = []
    for path in pool:
        name = path.name.lower()
        if patterns and not any(fnmatch.fnmatch(name, pattern) for pattern in patterns):
            continue
        try:
            header_df = read_table_file(path, nrows=0)
        except Exception as exc:  # pragma: no cover - defensive path
            warnings.append(f"{logical_table}: could not inspect {path}: {exc}")
            continue
        _, mapping = standardize_columns(header_df, table_config)
        if has_required_columns(mapping, table_config):
            matched.append(path)
        else:
            missing = [
                col for col in table_config.get("required_columns", []) if col not in mapping
            ]
            warnings.append(
                f"{logical_table}: skipped {path.name}; missing required mapped columns {missing}"
            )
    return matched, warnings


def load_logical_table(
    raw_dir: Path, logical_table: str, table_config: dict[str, Any], use_demo_if_no_raw: bool
) -> tuple[pd.DataFrame, list[dict[str, Any]], list[str]]:
    files, warnings = detect_source_files(raw_dir, logical_table, table_config, use_demo_if_no_raw)
    frames: list[pd.DataFrame] = []
    source_records: list[dict[str, Any]] = []
    for path in files:
        try:
            raw = read_table_file(path)
            mapped, mapping = standardize_columns(raw, table_config)
        except Exception as exc:
            warnings.append(f"{logical_table}: failed to load {path}: {exc}")
            continue
        if mapped.empty and len(raw) > 0:
            warnings.append(f"{logical_table}: no configured columns mapped from {path.name}")
            continue
        frames.append(mapped)
        source_records.append(
            {
                "logical_table": logical_table,
                "source_file": str(path),
                "rows_loaded": int(len(mapped)),
                "column_mapping": json.dumps(mapping, sort_keys=True),
            }
        )
    if not frames:
        return pd.DataFrame(), source_records, warnings
    combined = pd.concat(frames, ignore_index=True, sort=False)
    return combined, source_records, warnings
