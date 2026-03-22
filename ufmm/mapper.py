from dataclasses import asdict
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

from .config import UFMMConfig
from .data import SequenceRecord, curate_records
from .motifs import classify_function, discover_candidates


def _bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        running = min(running, val)
        q[i] = running
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def _group_persistence(group_hits: Dict[str, int], min_group_prevalence: float) -> float:
    values = np.array(list(group_hits.values()), dtype=float)
    if values.size == 0:
        return 0.0
    nonzero = np.mean(values > 0)
    abundance = np.mean(values / (values.max() + 1e-9))
    gate = float(nonzero >= min_group_prevalence)
    return float(gate * (0.65 * nonzero + 0.35 * abundance))


def _permutation_pvalue(observed: float, values: np.ndarray, rounds: int, rng: np.random.Generator) -> float:
    if values.size <= 1:
        return 1.0
    hits = 1
    for _ in range(rounds):
        perm = rng.permutation(values)
        perm_score = np.mean(perm > 0)
        if perm_score >= observed:
            hits += 1
    return hits / (rounds + 1)


def build_motif_function_map(records: Iterable[SequenceRecord], config: UFMMConfig) -> pd.DataFrame:
    curated = curate_records(records, config.min_length, config.max_length)
    if not curated:
        return pd.DataFrame()

    candidates = discover_candidates(
        records=curated,
        kmer_min=config.kmer_min,
        kmer_max=config.kmer_max,
        min_global_occurrences=config.min_global_occurrences,
        max_motifs=config.max_motifs,
    )

    rng = np.random.default_rng(config.random_seed)
    rows = []
    for motif, count, group_hits in candidates:
        vals = np.array(list(group_hits.values()), dtype=float)
        persistence = _group_persistence(group_hits, config.min_group_prevalence)
        pval = _permutation_pvalue(float(np.mean(vals > 0)), vals, config.permutation_rounds, rng)
        function = classify_function(motif)
        prediction = (
            f"Mutational scan of motif '{motif}' in proteins annotated as {function}; "
            "measure replication, regulation, and host-interaction phenotypes."
        )
        rows.append(
            {
                "motif": motif,
                "global_occurrences": int(count),
                "groups_supported": int(np.sum(vals > 0)),
                "group_hits": group_hits,
                "functional_class": function,
                "functional_conservation_score": round(persistence, 4),
                "p_value": pval,
                "prediction": prediction,
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df["q_value"] = _bh_qvalues(df["p_value"].to_numpy())
    df = df.sort_values(["q_value", "functional_conservation_score", "global_occurrences"], ascending=[True, False, False])
    df["passes_validation"] = df["q_value"] <= config.qvalue_threshold
    df["config"] = str(asdict(config))
    return df.reset_index(drop=True)
