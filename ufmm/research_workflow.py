"""Research-grade orchestration utilities for UFMM.

This module adds a reproducible end-to-end workflow that can be used directly in
papers, supplementary methods, and collaborative datasets.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
import requests

from .config import UFMMConfig
from .data import SequenceRecord
from .mapper import build_motif_function_map

UNIPROT_FASTA_ENDPOINT = "https://rest.uniprot.org/uniprotkb/stream"


@dataclass(frozen=True)
class RetrievalSpec:
    """Describes a single biological source query for reproducible retrieval."""

    group: str
    taxa: tuple[str, ...]
    limit_per_taxon: int = 100
    reviewed_only: bool = True
    uniprot_keyword: str | None = "KW-1194"


@dataclass(frozen=True)
class WorkflowSummary:
    """Machine-readable summary metadata produced by :func:`run_research_workflow`."""

    generated_at_utc: str
    ufmm_config: dict
    retrieval_specs: list[dict]
    total_sequences_loaded: int
    motifs_discovered: int
    motifs_validated: int
    output_csv: str


def _build_uniprot_query(spec: RetrievalSpec, taxon: str) -> str:
    clauses = [f"organism:{taxon}"]
    if spec.reviewed_only:
        clauses.append("reviewed:true")
    if spec.uniprot_keyword:
        clauses.append(f"keyword:{spec.uniprot_keyword}")
    return " AND ".join(clauses)


def fetch_uniprot_records(spec: RetrievalSpec, timeout_seconds: int = 30) -> list[SequenceRecord]:
    """Fetch amino-acid sequences from UniProt for a retrieval specification."""

    rows: list[SequenceRecord] = []
    for taxon in spec.taxa:
        query = _build_uniprot_query(spec, taxon)
        response = requests.get(
            UNIPROT_FASTA_ENDPOINT,
            params={"format": "fasta", "query": query, "size": spec.limit_per_taxon},
            timeout=timeout_seconds,
        )
        response.raise_for_status()

        header = None
        seq_parts: list[str] = []
        for line in response.text.splitlines():
            if line.startswith(">"):
                if header and seq_parts:
                    rows.append(
                        SequenceRecord(
                            seq_id=header.split()[0],
                            sequence="".join(seq_parts),
                            group=spec.group,
                            source_file=f"uniprot:{taxon}",
                        )
                    )
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

        if header and seq_parts:
            rows.append(
                SequenceRecord(
                    seq_id=header.split()[0],
                    sequence="".join(seq_parts),
                    group=spec.group,
                    source_file=f"uniprot:{taxon}",
                )
            )

    return rows


def save_research_outputs(
    dataframe: pd.DataFrame,
    summary: WorkflowSummary,
    out_csv: Path,
    out_summary_json: Path,
    out_methods_md: Path,
) -> None:
    """Save UFMM outputs as publication-friendly artifacts."""

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_summary_json.parent.mkdir(parents=True, exist_ok=True)
    out_methods_md.parent.mkdir(parents=True, exist_ok=True)

    dataframe.to_csv(out_csv, index=False)
    out_summary_json.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")

    methods = [
        "# UFMM Computational Methods Snapshot",
        "",
        f"Generated: {summary.generated_at_utc}",
        f"Total curated sequences: {summary.total_sequences_loaded}",
        f"Candidate motifs discovered: {summary.motifs_discovered}",
        f"Validated motifs (q-value filter): {summary.motifs_validated}",
        "",
        "## Retrieval specifications",
    ]
    for spec in summary.retrieval_specs:
        methods.append(
            f"- {spec['group']}: taxa={', '.join(spec['taxa'])}; "
            f"limit_per_taxon={spec['limit_per_taxon']}; reviewed_only={spec['reviewed_only']}; "
            f"uniprot_keyword={spec['uniprot_keyword']}"
        )

    methods.extend(
        [
            "",
            "## UFMM configuration",
            "```json",
            json.dumps(summary.ufmm_config, indent=2),
            "```",
        ]
    )

    out_methods_md.write_text("\n".join(methods) + "\n", encoding="utf-8")


def run_research_workflow(
    retrieval_specs: Sequence[RetrievalSpec],
    config: UFMMConfig,
    out_csv: str | Path,
    out_summary_json: str | Path,
    out_methods_md: str | Path,
    additional_records: Iterable[SequenceRecord] | None = None,
) -> pd.DataFrame:
    """Execute a reproducible UFMM workflow from public retrieval specs."""

    all_records: list[SequenceRecord] = []
    for spec in retrieval_specs:
        all_records.extend(fetch_uniprot_records(spec))

    if additional_records is not None:
        all_records.extend(additional_records)

    df = build_motif_function_map(all_records, config=config)

    summary = WorkflowSummary(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        ufmm_config=asdict(config),
        retrieval_specs=[asdict(spec) for spec in retrieval_specs],
        total_sequences_loaded=len(all_records),
        motifs_discovered=len(df),
        motifs_validated=int(df.get("passes_validation", pd.Series(dtype=bool)).sum()) if not df.empty else 0,
        output_csv=str(Path(out_csv)),
    )

    save_research_outputs(
        dataframe=df,
        summary=summary,
        out_csv=Path(out_csv),
        out_summary_json=Path(out_summary_json),
        out_methods_md=Path(out_methods_md),
    )
    return df
