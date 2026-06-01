#!/usr/bin/env python3
"""Validate drug mention extraction output from the oncology R pipeline."""

from pathlib import Path
import csv
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "05_drug_therapy_extraction" / "oncology_drug_mentions.csv"
PIPELINE = ROOT / "13_R_scripts" / "01_oncology_analysis_pipeline.R"


def count_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> int:
    if not OUT.exists():
        result = subprocess.run(["Rscript", str(PIPELINE)], check=False)
        if result.returncode != 0:
            return result.returncode
    if not OUT.exists():
        print("Drug mention output was not generated.", file=sys.stderr)
        return 1
    print(f"drug_mentions_rows={count_rows(OUT)}")
    print(f"output={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

