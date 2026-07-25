# botas/core/scoring.py
from __future__ import annotations

import math


def score_to_mapq(score: int, *, read_len: int, num_windows: int) -> int:
    """
    Simple MAPQ heuristic for early development.

    - Higher parasail score → higher MAPQ
    - More candidate windows → lower MAPQ (more ambiguity)

    This is deliberately conservative and easy to replace later with
    calibration (e.g., score gap between best/second-best).
    """
    if score <= 0 or read_len <= 0:
        return 0

    # Normalize score by read length
    s = score / max(1.0, float(read_len))

    # Penalize ambiguity
    amb_pen = math.log2(max(2, num_windows))  # 1 window ~ low penalty, many windows ~ higher penalty

    # Convert to 0..60 scale (rough)
    mq = int(max(0.0, min(60.0, (s * 20.0) - (amb_pen * 5.0))))
    return mq
