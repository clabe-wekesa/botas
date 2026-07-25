# botas/core/cigar.py
"""
CIGAR utilities (single source of truth for ref/query consumption helpers).

We keep these here so alignment modules (align_core / align_pe) and evaluators
do not re-implement CIGAR parsing logic.
"""

from __future__ import annotations

import re
from typing import Iterable, List, Tuple

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _iter_cigar_ops(cigar: str) -> Iterable[Tuple[int, str]]:
    """
    Yield (length, op) for each CIGAR operation.
    Accepts standard ops: MIDNSHP=X.
    """
    if not cigar:
        return
    for n, op in _CIGAR_RE.findall(cigar):
        yield int(n), op


def ref_aligned_length(cigar: str) -> int:
    """
    Reference-aligned length consumed by alignment operations on the reference.
    Counts ops that consume REF: M, D, N, =, X
    Does NOT count I, S, H, P.
    """
    ref_len = 0
    for n, op in _iter_cigar_ops(cigar):
        if op in ("M", "D", "N", "=", "X"):
            ref_len += n
    return ref_len


def ref_leading_consumption(cigar: str) -> int:
    """
    Leading reference consumption *before* the first ref-consuming aligned op.

    This is useful when aligners/padders introduce leading operations and you need
    to know how much REF would be shifted by the leading part.

    Rules:
      - Leading ops that consume REF: D, N contribute.
      - Leading ops that do NOT consume REF: S, H, I, P contribute 0.
      - Stop at first aligned ref-consuming op: M, =, X (and also D/N if they appear
        after alignment begins, but for "leading" we stop at the first M/= /X).
    """
    lead = 0
    for n, op in _iter_cigar_ops(cigar):
        if op in ("M", "=", "X"):
            break
        if op in ("D", "N"):
            lead += n
        # I/S/H/P: do not consume REF, ignore
    return lead


def trim_ref_padding_from_cigar(
    cigar: str,
    left_ref_pad: int,
    right_ref_pad: int,
) -> str:
    """
    Trim *reference padding* that was introduced by slicing the reference with pads.

    Scenario (common in PE rescue):
      - You slice reference around a center with pad bases on both sides.
      - You align a read to that slice.
      - The resulting CIGAR is in the slice coordinate system.
      - To map back cleanly (and avoid "aligning into padding"), you can trim away
        ref-only consumption that falls in the artificial pads.

    We trim from the CIGAR ends by removing REF-consuming ops that lie fully/partly
    in the pads. This function operates on the REF axis, not the read axis.

    Notes:
      - We only trim against REF-consuming ops: M, D, N, =, X.
      - If trimming splits an op, we shorten it.
      - Ops that do not consume REF (I/S/H/P) are preserved unless they become
        stranded at ends after full trim; we keep them as-is to avoid surprises.

    Returns a new CIGAR string (may be empty if everything trimmed).
    """
    if not cigar or (left_ref_pad <= 0 and right_ref_pad <= 0):
        return cigar

    ops: List[Tuple[int, str]] = list(_iter_cigar_ops(cigar))
    if not ops:
        return cigar

    def consumes_ref(op: str) -> bool:
        return op in ("M", "D", "N", "=", "X")

    # Trim left pad
    ltrim = max(0, left_ref_pad)
    i = 0
    while ltrim > 0 and i < len(ops):
        n, op = ops[i]
        if not consumes_ref(op):
            i += 1
            continue
        if n <= ltrim:
            ltrim -= n
            ops[i] = (0, op)
            i += 1
        else:
            ops[i] = (n - ltrim, op)
            ltrim = 0

    # Trim right pad
    rtrim = max(0, right_ref_pad)
    j = len(ops) - 1
    while rtrim > 0 and j >= 0:
        n, op = ops[j]
        if n == 0:
            j -= 1
            continue
        if not consumes_ref(op):
            j -= 1
            continue
        if n <= rtrim:
            rtrim -= n
            ops[j] = (0, op)
            j -= 1
        else:
            ops[j] = (n - rtrim, op)
            rtrim = 0

    # Rebuild (drop zero-len ops)
    out = []
    for n, op in ops:
        if n > 0:
            out.append(f"{n}{op}")
    return "".join(out)
