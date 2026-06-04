from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_database import build_database


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate raw atlas inputs without writing DuckDB output.")
    parser.add_argument("--project-dir", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--raw-dir", type=Path, default=None)
    args = parser.parse_args()
    result = build_database(
        project_dir=args.project_dir,
        raw_dir=args.raw_dir,
        write_db=False,
        export_outputs=False,
    )
    print(f"Validation report: {result['validation_report']}")
    for warning in result["warnings"]:
        print(f"{warning['severity']}: {warning['message']}")


if __name__ == "__main__":
    main()
