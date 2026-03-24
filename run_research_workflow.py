#!/usr/bin/env python3
"""Publication-oriented runner for the UFMM research workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from ufmm import UFMMConfig
from ufmm.research_workflow import RetrievalSpec, run_research_workflow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a reproducible UniProt-to-UFMM workflow")
    parser.add_argument("--bacteria", nargs="*", default=[])
    parser.add_argument("--viruses", nargs="*", default=[])
    parser.add_argument("--mges", nargs="*", default=[])
    parser.add_argument("--host", nargs="*", default=[])

    parser.add_argument("--limit-per-taxon", type=int, default=100)
    parser.add_argument("--keyword", type=str, default="KW-1194")

    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    parser.add_argument("--output-methods", type=Path, required=True)

    parser.add_argument("--kmer-min", type=int, default=6)
    parser.add_argument("--kmer-max", type=int, default=12)
    parser.add_argument("--min-global-occurrences", type=int, default=8)
    parser.add_argument("--permutation-rounds", type=int, default=200)
    parser.add_argument("--qvalue-threshold", type=float, default=0.1)
    return parser.parse_args()


def _add_group_specs(storage: list[RetrievalSpec], group: str, taxa: list[str], limit: int, keyword: str) -> None:
    cleaned = tuple(t.strip() for t in taxa if t and t.strip())
    if cleaned:
        storage.append(
            RetrievalSpec(
                group=group,
                taxa=cleaned,
                limit_per_taxon=limit,
                uniprot_keyword=keyword,
            )
        )


def main() -> None:
    args = parse_args()
    specs: list[RetrievalSpec] = []
    _add_group_specs(specs, "Bacteria", args.bacteria, args.limit_per_taxon, args.keyword)
    _add_group_specs(specs, "Virus", args.viruses, args.limit_per_taxon, args.keyword)
    _add_group_specs(specs, "MGE", args.mges, args.limit_per_taxon, args.keyword)
    _add_group_specs(specs, "Host", args.host, args.limit_per_taxon, args.keyword)

    if not specs:
        raise SystemExit("No taxa were provided. Use --bacteria/--viruses/--mges/--host.")

    config = UFMMConfig(
        kmer_min=args.kmer_min,
        kmer_max=args.kmer_max,
        min_global_occurrences=args.min_global_occurrences,
        permutation_rounds=args.permutation_rounds,
        qvalue_threshold=args.qvalue_threshold,
    )

    df = run_research_workflow(
        retrieval_specs=specs,
        config=config,
        out_csv=args.output_csv,
        out_summary_json=args.output_summary,
        out_methods_md=args.output_methods,
    )

    print(f"Wrote {len(df)} motifs to {args.output_csv}")
    if not df.empty:
        print(f"Validated motifs: {int(df['passes_validation'].sum())}")


if __name__ == "__main__":
    main()
