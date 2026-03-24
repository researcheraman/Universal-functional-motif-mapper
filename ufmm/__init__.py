"""Universal Functional Motif Mapper (UFMM) core package."""

from .config import UFMMConfig
from .data import SequenceRecord, load_fasta_group, curate_records
from .mapper import build_motif_function_map
from .research_workflow import RetrievalSpec, WorkflowSummary, run_research_workflow

__all__ = [
    "UFMMConfig",
    "SequenceRecord",
    "load_fasta_group",
    "curate_records",
    "build_motif_function_map",
    "RetrievalSpec",
    "WorkflowSummary",
    "run_research_workflow",
]
