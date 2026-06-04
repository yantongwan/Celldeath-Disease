from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.normalization import normalize_death_mode, normalize_disease, normalize_drug_name, normalize_gene_symbol, normalize_pmid


NORMALIZERS = {
    "pmid": normalize_pmid,
    "death_mode": normalize_death_mode,
    "disease": normalize_disease,
    "gene": normalize_gene_symbol,
    "drug": normalize_drug_name,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize a single CSV/TSV column for audit use.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--column", required=True)
    parser.add_argument("--kind", choices=sorted(NORMALIZERS), required=True)
    args = parser.parse_args()

    sep = "\t" if args.input.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(args.input, sep=sep, dtype=str, keep_default_na=False)
    if args.column not in df.columns:
        raise SystemExit(f"Column not found: {args.column}")
    df[f"{args.column}_normalized"] = df[args.column].map(NORMALIZERS[args.kind])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
