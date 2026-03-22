#!/usr/bin/env python3
"""CLI runner for Universal Functional Motif Mapper.

Example:
python run_ufmm.py \
  --bacteria data/bacteria.faa \
  --viruses data/viruses.faa \
  --mges data/mges.faa \
  --output results/motif_function_map.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ufmm import UFMMConfig, build_motif_function_map, load_fasta_group


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Universal Functional Motif Mapper")
    p.add_argument("--bacteria", type=Path, nargs="*", default=[])
    p.add_argument("--viruses", type=Path, nargs="*", default=[])
    p.add_argument("--mges", type=Path, nargs="*", default=[])
    p.add_argument("--host", type=Path, nargs="*", default=[])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--kmer-min", type=int, default=6)
    p.add_argument("--kmer-max", type=int, default=12)
    p.add_argument("--min-global-occurrences", type=int, default=8)
    p.add_argument("--permutation-rounds", type=int, default=200)
    p.add_argument("--qvalue-threshold", type=float, default=0.1)
    return p.parse_args()


def _load(paths: list[Path], group: str):
    recs = []
    for path in paths:
        recs.extend(load_fasta_group(path, group=group))
    return recs


def main() -> None:
    args = parse_args()

    records = []
    records.extend(_load(args.bacteria, "Bacteria"))
    records.extend(_load(args.viruses, "Virus"))
    records.extend(_load(args.mges, "MGE"))
    records.extend(_load(args.host, "Host"))

    config = UFMMConfig(
        kmer_min=args.kmer_min,
        kmer_max=args.kmer_max,
        min_global_occurrences=args.min_global_occurrences,
        permutation_rounds=args.permutation_rounds,
        qvalue_threshold=args.qvalue_threshold,
    )

    df = build_motif_function_map(records, config=config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    print(f"Wrote {len(df)} motifs to {args.output}")
    if not df.empty:
        validated = int(df["passes_validation"].sum())
        print(f"Validated motifs (q <= {config.qvalue_threshold}): {validated}")


if __name__ == "__main__":
    main()
