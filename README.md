# Universal Functional Motif Mapper (UFMM)

This repository now includes a reusable computational pipeline to discover and validate
function-first motifs across bacteria, viruses, mobile genetic elements (MGEs), and host-linked datasets.

## What it does

- Loads FASTA datasets grouped by biological system.
- Curates sequences (alphabet cleanup + length filters).
- Discovers recurrent motifs across a configurable k-mer range.
- Scores **functional conservation** by cross-group persistence.
- Estimates empirical p-values with permutation testing.
- Controls false discoveries with Benjamini-Hochberg q-values.
- Emits a motif-function map with testable predictions.

## Quickstart

```bash
python run_ufmm.py \
  --bacteria data/bacteria.faa \
  --viruses data/viruses.faa \
  --mges data/mges.faa \
  --output results/motif_function_map.csv
```

## Output columns

- `motif`
- `global_occurrences`
- `groups_supported`
- `group_hits`
- `functional_class`
- `functional_conservation_score`
- `p_value`
- `q_value`
- `passes_validation`
- `prediction`
- `config`

## Notes

This is a rigorous computational starting point intended to generate hypotheses for wet-lab validation
(e.g., mutational scans, motif swaps, replication/regulation assays).
