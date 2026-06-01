#!/usr/bin/env python3
"""Comorbidity PubMed pipeline

This script merges two stages:
1) Filter disease terms by PubMed hit count (precision query).
2) Build disease–disease co-mention network (all pairs totals) and, for high-count pairs,
   extract yearly trends.

It is a refactor of the user's original scripts:
- 0_filter_by_count.py
- 0_pubmed_run_all.py

Key improvements:
- No hard-coded Entrez credentials; use CLI args or environment variables.
- Safer retries/backoff for transient NCBI/HTTP errors.
- Streaming combinations (no huge list allocation).
- Consistent query builder and common utilities.

Usage examples:
  python comorbidity_pubmed_pipeline.py filter --input diseases.txt --min-count 10000
  python comorbidity_pubmed_pipeline.py network --disease-list diseases_top.txt --trend-threshold 500 \
      --start-year 2014 --end-year 2024 --start-index 0

Environment variables:
  ENTREZ_EMAIL, ENTREZ_API_KEY
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
import ssl
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from Bio import Entrez
from urllib.error import HTTPError, URLError


# ---------------------------
# SSL workaround (kept from original)
# ---------------------------
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    _create_unverified_https_context = None
else:
    ssl._create_default_https_context = _create_unverified_https_context


# ---------------------------
# Utilities
# ---------------------------

def configure_entrez(email: str, api_key: Optional[str] = None) -> None:
    """Configure NCBI Entrez. Email is required by NCBI policy."""
    if not email:
        raise ValueError("Entrez email is required. Pass --email or set ENTREZ_EMAIL.")
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def build_precision_query(term: str) -> str:
    """Precision query used in both stages.

    Note:
      - We keep the user's intent: include both controlled vocabulary (MeSH Terms) and
        free text in title/abstract.
      - PubMed field tag equivalent for Title/Abstract is [tiab].
    """
    term = term.strip()
    # Prefer tiab (compact) but keep Title/Abstract spelling acceptable too.
    return f'"{term}"[MeSH Terms] OR "{term}"[tiab]'


def entrez_esearch_count(term: str, db: str = "pubmed", max_tries: int = 6, base_sleep: float = 0.4) -> int:
    """Return PubMed count for a query, with retries/backoff.

    Returns:
      count (>=0) on success, -1 on failure.
    """
    # Exponential backoff with jitter
    for attempt in range(1, max_tries + 1):
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=0)
            record = Entrez.read(handle)
            handle.close()
            return int(record.get("Count", 0))
        except HTTPError as e:
            # NCBI can return 429/5xx; backoff then retry.
            code = getattr(e, "code", None)
            sleep_s = base_sleep * (2 ** (attempt - 1))
            # Slight extra delay for rate limiting
            if code == 429:
                sleep_s *= 2
            print(f" [HTTPError {code} try {attempt}/{max_tries}]", end="", flush=True)
            time.sleep(sleep_s)
        except (URLError, OSError, RuntimeError, ValueError) as e:
            sleep_s = base_sleep * (2 ** (attempt - 1))
            print(f" [Error try {attempt}/{max_tries}: {e}]", end="", flush=True)
            time.sleep(sleep_s)
        except Exception as e:
            # Unknown: still retry a couple times.
            sleep_s = base_sleep * (2 ** (attempt - 1))
            print(f" [Unexpected try {attempt}/{max_tries}: {e}]", end="", flush=True)
            time.sleep(sleep_s)
    return -1


def read_terms_from_txt(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        terms = [line.strip() for line in f if line.strip()]
    # de-dup and stable ordering
    return sorted(set(terms))


def read_terms_from_csv_first_col(path: str, skip_header: bool = True) -> List[str]:
    terms: List[str] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        if skip_header:
            next(reader, None)
        for row in reader:
            if not row:
                continue
            t = row[0].strip()
            if t:
                terms.append(t)
    return sorted(set(terms))


def iter_pairs_stream(diseases: Sequence[str]) -> Iterator[Tuple[str, str]]:
    # stream, do not materialize list
    return itertools.combinations(diseases, 2)


def ensure_parent_dir(path: str) -> None:
    d = os.path.dirname(os.path.abspath(path))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


# ---------------------------
# Stage 1: filter
# ---------------------------

def run_filter(
    input_file: str,
    output_list: str,
    output_csv: str,
    min_paper_count: int,
    sleep_s: float,
) -> None:
    diseases = read_terms_from_txt(input_file)
    total = len(diseases)

    print("=== Stage 1: Filter disease terms by PubMed count ===")
    print(f"Input: {input_file}")
    print(f"Terms: {total}")
    print(f"Precision query: MeSH Terms OR tiab")
    print(f"Threshold: >= {min_paper_count}")
    print("-" * 60)

    kept: List[str] = []
    stats: List[Tuple[str, int]] = []

    for i, disease in enumerate(diseases, start=1):
        print(f"[{i}/{total}] {disease:<45} ... ", end="", flush=True)
        q = build_precision_query(disease)
        c = entrez_esearch_count(q)
        if c >= min_paper_count:
            print(f"✅ keep ({c})")
            kept.append(disease)
            stats.append((disease, c))
        else:
            if c == -1:
                print("❌ error")
            else:
                print(f"❌ drop ({c})")
        time.sleep(sleep_s)

    ensure_parent_dir(output_list)
    ensure_parent_dir(output_csv)

    with open(output_list, "w", encoding="utf-8") as f:
        f.write("\n".join(kept))

    with open(output_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Disease", "Article_Count"])
        w.writerows(stats)

    print("\n" + "=" * 60)
    print("Stage 1 done")
    print(f"Original: {total}")
    print(f"Kept: {len(kept)}")
    print(f"Saved list: {output_list}")
    print(f"Saved counts: {output_csv}")


# ---------------------------
# Stage 2: network
# ---------------------------

def load_diseases_auto(disease_list: Optional[str], disease_csv: Optional[str]) -> List[str]:
    if disease_list and os.path.exists(disease_list):
        print(f"📖 Loading terms from {disease_list}")
        return read_terms_from_txt(disease_list)
    if disease_csv and os.path.exists(disease_csv):
        print(f"📖 Loading terms from {disease_csv} (first column)")
        return read_terms_from_csv_first_col(disease_csv)
    raise FileNotFoundError(
        "No disease input found. Provide --disease-list (txt) or --disease-csv (csv)."
    )


def yearly_counts_for_pair(full_query: str, start_year: int, end_year: int, sleep_s: float) -> Optional[Dict[int, int]]:
    yearly: Dict[int, int] = {}
    for y in range(start_year, end_year + 1):
        # PubMed supports year-only PDAT; keep simple and fast
        qy = f"{full_query} AND {y}[pdat]"
        c = entrez_esearch_count(qy)
        if c == -1:
            return None
        yearly[y] = c
        time.sleep(sleep_s)
    return yearly


def run_network(
    disease_list: Optional[str],
    disease_csv: Optional[str],
    output_totals: str,
    output_trends: str,
    trend_threshold: int,
    start_year: int,
    end_year: int,
    start_index: int,
    sleep_s: float,
) -> None:
    diseases = load_diseases_auto(disease_list, disease_csv)
    n = len(diseases)
    if n < 2:
        print("Need at least 2 diseases.")
        return

    # total pairs count without materializing the list
    total_pairs = n * (n - 1) // 2

    print("=== Stage 2: Build disease–disease PubMed co-mention network ===")
    print(f"Diseases: {n}")
    print(f"Pairs: {total_pairs}")
    print(f"Totals output: {output_totals}")
    print(f"Trends output: {output_trends} (only totals > {trend_threshold})")
    print(f"Years: {start_year}–{end_year}")
    print(f"Resume from pair index: {start_index}")
    print("-" * 60)

    ensure_parent_dir(output_totals)
    ensure_parent_dir(output_trends)

    mode = "w" if start_index == 0 else "a"

    with open(output_totals, mode, encoding="utf-8", newline="") as f_total, open(
        output_trends, mode, encoding="utf-8", newline=""
    ) as f_trend:
        w_total = csv.writer(f_total)
        w_trend = csv.writer(f_trend)

        if start_index == 0:
            w_total.writerow(["Disease_A", "Disease_B", "Total_Count"])
            header = ["Disease_A", "Disease_B", "Total_Count"] + [str(y) for y in range(start_year, end_year + 1)]
            w_trend.writerow(header)

        # streaming enumeration
        for idx, (d1, d2) in enumerate(iter_pairs_stream(diseases)):
            if idx < start_index:
                continue

            print(f"[{idx+1}/{total_pairs}] {d1[:22]}.. <-> {d2[:22]}.. ", end="", flush=True)

            q1 = build_precision_query(d1)
            q2 = build_precision_query(d2)
            full = f"({q1}) AND ({q2})"

            total = entrez_esearch_count(full)
            if total == -1:
                print("❌ total error")
                continue

            w_total.writerow([d1, d2, total])
            f_total.flush()

            if total > trend_threshold:
                print(f"🔥 total={total} extracting trends... ", end="", flush=True)
                trends = yearly_counts_for_pair(full, start_year, end_year, sleep_s)
                if trends is None:
                    print("❌ trend error")
                else:
                    row = [d1, d2, total] + [trends[y] for y in range(start_year, end_year + 1)]
                    w_trend.writerow(row)
                    f_trend.flush()
                    print("✅")
            else:
                print(f"⚪ total={total} skip")

            time.sleep(sleep_s)

    print("\n" + "=" * 60)
    print("Stage 2 done")
    print(f"Totals: {output_totals}")
    print(f"Trends: {output_trends}")


# ---------------------------
# CLI
# ---------------------------

def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="PubMed comorbidity pipeline: filter disease list and build co-mention network."
    )
    p.add_argument("--email", default=os.getenv("ENTREZ_EMAIL", ""), help="NCBI Entrez email (or set ENTREZ_EMAIL)")
    p.add_argument("--api-key", default=os.getenv("ENTREZ_API_KEY"), help="NCBI API key (or set ENTREZ_API_KEY)")
    p.add_argument("--sleep", type=float, default=0.12, help="Base sleep between requests (seconds)")

    sub = p.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("filter", help="Filter disease terms by PubMed count")
    p1.add_argument("--input", default="diseases.txt", help="Input txt list (one term per line)")
    p1.add_argument("--output-list", default="diseases_top.txt", help="Output txt (kept terms)")
    p1.add_argument("--output-csv", default="diseases_counts.csv", help="Output csv (term,count)")
    p1.add_argument("--min-count", type=int, default=10000, help="Minimum PubMed count to keep")

    p2 = sub.add_parser("network", help="Build disease pair co-mention network")
    p2.add_argument("--disease-list", default="diseases_top.txt", help="Input txt list (preferred if exists)")
    p2.add_argument("--disease-csv", default="diseases_counts.csv", help="Fallback csv if list not present")
    p2.add_argument("--output-totals", default="disease_network_all_totals.csv")
    p2.add_argument("--output-trends", default="disease_network_trends_gt500.csv")
    p2.add_argument("--trend-threshold", type=int, default=500)
    p2.add_argument("--start-year", type=int, default=2014)
    p2.add_argument("--end-year", type=int, default=2024)
    p2.add_argument("--start-index", type=int, default=0, help="Resume index (0-based pair index)")

    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_parser().parse_args(argv)

    try:
        configure_entrez(args.email, args.api_key)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 2

    if args.cmd == "filter":
        if not os.path.exists(args.input):
            print(f"Input file not found: {args.input}", file=sys.stderr)
            return 2
        run_filter(
            input_file=args.input,
            output_list=args.output_list,
            output_csv=args.output_csv,
            min_paper_count=args.min_count,
            sleep_s=args.sleep,
        )
        return 0

    if args.cmd == "network":
        run_network(
            disease_list=args.disease_list,
            disease_csv=args.disease_csv,
            output_totals=args.output_totals,
            output_trends=args.output_trends,
            trend_threshold=args.trend_threshold,
            start_year=args.start_year,
            end_year=args.end_year,
            start_index=args.start_index,
            sleep_s=args.sleep,
        )
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
