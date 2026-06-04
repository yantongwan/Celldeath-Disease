from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_database import export_static_tables


def main() -> None:
    parser = argparse.ArgumentParser(description="Export DuckDB tables and views to CSV and Parquet.")
    parser.add_argument("--db-path", type=Path, default=PROJECT_ROOT / "data" / "processed" / "rcd_atlas.duckdb")
    parser.add_argument("--exports-dir", type=Path, default=PROJECT_ROOT / "data" / "exports")
    parser.add_argument("--parquet-dir", type=Path, default=PROJECT_ROOT / "data" / "processed" / "rcd_atlas.parquet")
    args = parser.parse_args()
    export_static_tables(args.db_path, args.exports_dir, args.parquet_dir)
    print(f"Exported CSV to {args.exports_dir}")
    print(f"Exported Parquet to {args.parquet_dir}")


if __name__ == "__main__":
    main()
