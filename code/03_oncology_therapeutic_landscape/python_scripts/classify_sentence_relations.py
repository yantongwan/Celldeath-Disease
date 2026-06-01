#!/usr/bin/env python3
"""Validate sentence-level relation output from the oncology R pipeline."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "09_sentence_level_relation_audit" / "drug_death_relation_sentences.csv"
SAMPLE = ROOT / "09_sentence_level_relation_audit" / "drug_relation_manual_validation_sample.csv"
PIPELINE = ROOT / "13_R_scripts" / "01_oncology_analysis_pipeline.R"


def count_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> int:
    if not OUT.exists() or not SAMPLE.exists():
        result = subprocess.run(["Rscript", str(PIPELINE)], check=False)
        if result.returncode != 0:
            return result.returncode
    if not OUT.exists():
        print("Sentence relation output was not generated.", file=sys.stderr)
        return 1
    print(f"sentence_relation_rows={count_rows(OUT)}")
    print(f"manual_validation_sample_rows={count_rows(SAMPLE) if SAMPLE.exists() else 0}")
    print(f"output={OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

