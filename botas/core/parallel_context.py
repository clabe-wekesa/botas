# botas/core/parallel_context.py
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class PEContext:
    """
    Paired-end execution context (multi-contig aware).

    Used by:
      - botas/core/pe_pool.py
      - botas/cli/align.py
    """
    ref_fa: str
    circular_all: bool
    circular_overhang_percent: int
    k: int
    step: int
    pad: int
    max_windows: int
    min_seed_hits: int
    max_insert: int
    expected_insert: int
    rescue_pad: int
    do_rescue: bool
    sensitive: bool
    debug_pairs: int
    prebuilt_contigs: object | None = None