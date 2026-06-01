#!/usr/bin/env python3
"""File-level validation for the oncology therapeutic death landscape module."""

from pathlib import Path
import csv
import sys


ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "02_oncology_scope_and_tumor_mapping/oncology_tumor_term_mapping.csv",
    "03_oncology_death_mode_hotness/oncology_death_mode_hotness.csv",
    "04_tumor_specific_death_profiles/tumor_death_selectivity.csv",
    "05_drug_therapy_extraction/drug_lexicon_normalized.csv",
    "05_drug_therapy_extraction/oncology_drug_mentions.csv",
    "09_sentence_level_relation_audit/drug_death_relation_sentences.csv",
    "09_sentence_level_relation_audit/drug_relation_manual_validation_sample.csv",
    "15_QC_reports/oncology_analysis_gate_decision_table.csv",
]


def count_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> int:
    rows = []
    ok = True
    for rel in REQUIRED:
        path = ROOT / rel
        exists = path.exists()
        ok = ok and exists
        rows.append({
            "file": rel,
            "exists_yes_no": "yes" if exists else "no",
            "row_count": count_rows(path) if exists else "",
        })
    out = ROOT / "15_QC_reports" / "oncology_output_file_validation.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "exists_yes_no", "row_count"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"validation_report={out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

