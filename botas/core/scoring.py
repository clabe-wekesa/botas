# botas/core/scoring.py
from __future__ import annotations

from typing import Optional


def score_to_mapq(
    score: int,
    *,
    read_len: int,
    second_score: Optional[int] = None,
) -> int:
    """
    Convert BOTAS Edlib alignment scores to a conservative 0..60 MAPQ heuristic.

    BOTAS defines alignment score as ``score = -edit_distance``; therefore an
    exact match has score 0 and poorer alignments have increasingly negative
    scores.

    MAPQ combines:
      1. alignment quality, measured as 1 - edit_distance/read_len; and
      2. separation from the second-best *distinct mapping locus*.

    If a second-best locus has the same score as the best locus, MAPQ is 0.
    If no second distinct locus is available, separation is treated as maximal.

    This is a heuristic mapping-confidence score, not a calibrated Phred
    probability of mapping error.
    """
    if read_len <= 0:
        return 0

    best_edits = max(0, -int(score))
    edit_fraction = min(1.0, best_edits / float(read_len))
    alignment_quality = 1.0 - edit_fraction

    if alignment_quality <= 0.0:
        return 0

    if second_score is None:
        separation = 1.0
    else:
        gap = int(score) - int(second_score)
        if gap <= 0:
            return 0

        # A gap of 10% of read length (or more) receives full separation
        # credit; smaller gaps are scaled linearly.
        gap_scale = max(1.0, 0.10 * float(read_len))
        separation = min(1.0, gap / gap_scale)

    mq = round(60.0 * alignment_quality * separation)
    return int(max(0, min(60, mq)))
