from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO


AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class SequenceRecord:
    """Curated amino-acid sequence with provenance metadata."""

    seq_id: str
    sequence: str
    group: str
    source_file: str


def _clean_sequence(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if ch in AA_ALPHABET)


def load_fasta_group(path: str | Path, group: str) -> List[SequenceRecord]:
    records: list[SequenceRecord] = []
    p = Path(path)
    for rec in SeqIO.parse(str(p), "fasta"):
        records.append(
            SequenceRecord(
                seq_id=rec.id,
                sequence=str(rec.seq),
                group=group,
                source_file=str(p),
            )
        )
    return records


def curate_records(
    records: Iterable[SequenceRecord],
    min_length: int,
    max_length: int,
) -> List[SequenceRecord]:
    curated: list[SequenceRecord] = []
    for rec in records:
        cleaned = _clean_sequence(rec.sequence)
        if min_length <= len(cleaned) <= max_length:
            curated.append(
                SequenceRecord(
                    seq_id=rec.seq_id,
                    sequence=cleaned,
                    group=rec.group,
                    source_file=rec.source_file,
                )
            )
    return curated
