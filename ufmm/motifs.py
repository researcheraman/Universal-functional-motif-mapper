from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple

from .data import SequenceRecord


def _iter_kmers(seq: str, k: int):
    for i in range(0, len(seq) - k + 1):
        yield seq[i : i + k]


def collect_kmer_counts(records: Iterable[SequenceRecord], k: int) -> tuple[Counter, Dict[str, Counter]]:
    global_counts: Counter = Counter()
    by_group: dict[str, Counter] = defaultdict(Counter)

    for rec in records:
        seen_this_seq = set()
        for kmer in _iter_kmers(rec.sequence, k):
            global_counts[kmer] += 1
            if kmer not in seen_this_seq:
                by_group[rec.group][kmer] += 1
                seen_this_seq.add(kmer)
    return global_counts, by_group


def prevalence_by_group(by_group: Dict[str, Counter], kmer: str) -> dict[str, int]:
    return {group: counter.get(kmer, 0) for group, counter in by_group.items()}


def classify_function(kmer: str) -> str:
    """Simple biological prior; replace/extend with curated models as data grows."""
    if "GKT" in kmer or "PLOOP" in kmer:
        return "replication_or_nucleotide_binding"
    if sum(1 for aa in kmer if aa in "KR") / len(kmer) >= 0.35:
        return "regulation_or_nucleic_acid_interaction"
    if sum(1 for aa in kmer if aa in "DEH") / len(kmer) >= 0.35:
        return "catalytic_or_processing"
    if sum(1 for aa in kmer if aa in "AILVFWY") / len(kmer) >= 0.55:
        return "host_interaction_or_membrane_association"
    return "unresolved_core_function"


def discover_candidates(
    records: List[SequenceRecord],
    kmer_min: int,
    kmer_max: int,
    min_global_occurrences: int,
    max_motifs: int,
) -> List[Tuple[str, int, Dict[str, int]]]:
    all_candidates: list[tuple[str, int, dict[str, int]]] = []
    for k in range(kmer_min, kmer_max + 1):
        global_counts, by_group = collect_kmer_counts(records, k)
        for kmer, cnt in global_counts.items():
            if cnt >= min_global_occurrences:
                all_candidates.append((kmer, cnt, prevalence_by_group(by_group, kmer)))

    all_candidates.sort(key=lambda t: (t[1], len([v for v in t[2].values() if v > 0])), reverse=True)
    return all_candidates[:max_motifs]
