from dataclasses import dataclass


@dataclass(frozen=True)
class UFMMConfig:
    """Configurable parameters for motif discovery and validation."""

    kmer_min: int = 6
    kmer_max: int = 12
    min_length: int = 40
    max_length: int = 5000
    min_group_prevalence: float = 0.25
    min_global_occurrences: int = 8
    max_motifs: int = 200
    permutation_rounds: int = 200
    random_seed: int = 7
    qvalue_threshold: float = 0.1
