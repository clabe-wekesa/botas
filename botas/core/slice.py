# botas/core/slice.py
"""
Reference slicing utilities.

Single source of truth for extracting reference slices
for rescue / extension, supporting linear and circular genomes.
"""

from __future__ import annotations
from typing import Tuple

from botas.core.circular import norm_pos


def slice_reference(
    *,
    ref_seq: str,
    center0: int,
    read_len: int,
    pad: int,
    circular: bool,
) -> Tuple[str, int]:
    """
    Extract a reference slice centered around center0.

    Returns:
        (ref_slice, slice_start0)

    slice_start0 is the LINEAR coordinate corresponding to
    ref_slice[0] in the reference space.

    Notes:
      - For circular genomes, the reference may be conceptually
        doubled (ref + ref) to allow wrap-spanning slices.
      - All normalization is done here.
    """
    L = len(ref_seq)
    if L == 0:
        return "", 0

    # Normalize center position
    c0 = norm_pos(center0, L)

    span = read_len + pad

    if circular:
        # Double reference to allow wrap
        ref2 = ref_seq + ref_seq
        center = c0 + L

        start = max(0, center - pad)
        end = min(len(ref2), center + read_len + pad)

        if end <= start:
            return "", 0

        return ref2[start:end], start

    # Linear case
    start = max(0, c0 - pad)
    end = min(L, c0 + read_len + pad)

    if end <= start:
        return "", 0

    return ref_seq[start:end], start
