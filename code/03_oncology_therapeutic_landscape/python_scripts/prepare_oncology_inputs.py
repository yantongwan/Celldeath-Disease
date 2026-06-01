#!/usr/bin/env python3
"""Prepare oncology module inputs by invoking the reproducible R pipeline.

This wrapper intentionally does not draw figures. It exists for users who want a
Python entry point for input preparation while keeping all plotting in R.
"""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "13_R_scripts" / "01_oncology_analysis_pipeline.R"


def main() -> int:
    if not SCRIPT.exists():
        print(f"Missing R pipeline: {SCRIPT}", file=sys.stderr)
        return 1
    result = subprocess.run(["Rscript", str(SCRIPT)], check=False)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())

