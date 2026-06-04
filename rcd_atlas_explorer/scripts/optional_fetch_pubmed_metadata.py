from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Optional disabled-by-default PubMed metadata enrichment placeholder. "
            "The main app never requires network access."
        )
    )
    parser.add_argument("--pmid-file", type=Path, required=False)
    parser.add_argument("--output", type=Path, required=False)
    parser.add_argument(
        "--enabled",
        action="store_true",
        help="Acknowledge that this optional enrichment would require an internet connection.",
    )
    args = parser.parse_args()
    if not args.enabled:
        raise SystemExit(
            "PubMed metadata enrichment is disabled by default. "
            "Use --enabled and implement a local policy-compliant fetcher before running online."
        )
    raise SystemExit(
        "No online fetcher is bundled. Add project-approved PubMed retrieval code here if enrichment is needed."
    )


if __name__ == "__main__":
    main()
