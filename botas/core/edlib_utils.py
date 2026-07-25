# botas/core/edlib_utils.py
"""
Edlib utility wrappers (single source of truth).

All edlib access MUST go through this module.
No direct `import edlib` elsewhere.
"""

from __future__ import annotations
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)

try:
    import edlib  # type: ignore
    _EDLIB_AVAILABLE = True
    _EDLIB_IMPORT_ERROR: Optional[Exception] = None
except ImportError as e:  # pragma: no cover
    edlib = None
    _EDLIB_AVAILABLE = False
    _EDLIB_IMPORT_ERROR = e


def require_edlib() -> None:
    if not _EDLIB_AVAILABLE:  # pragma: no cover
        raise ImportError(
            "edlib is required. Install with: pip install edlib"
        ) from _EDLIB_IMPORT_ERROR


def edlib_align_hw_path(
    query: str,
    ref: str,
) -> Optional[Tuple[int, str, int, int]]:
    """
    Semi-global alignment (HW) with CIGAR + best location.

    Returns:
        (score, cigar, ref_beg, ref_end)

    FIX: Previously called edlib.align() TWICE — once with task="path" and
    once with task="locations". task="path" already returns locations in the
    result dict, so a single call is sufficient.
    """
    require_edlib()

    res = edlib.align(query, ref, mode="HW", task="path")
    if res["editDistance"] < 0 or not res.get("cigar") or not res.get("locations"):
        return None

    score = -int(res["editDistance"])
    cigar = str(res["cigar"])
    ref_beg, ref_end = res["locations"][0]

    return score, cigar, int(ref_beg), int(ref_end)


def edlib_align_hw_locations(
    query: str,
    ref: str,
) -> Optional[Tuple[int, int, int]]:
    """
    Fast HW alignment returning only (edits, ref_beg, ref_end).
    """
    require_edlib()

    res = edlib.align(query, ref, mode="HW", task="locations")
    if res["editDistance"] < 0 or not res.get("locations"):
        return None

    ref_beg, ref_end = res["locations"][0]
    return int(res["editDistance"]), int(ref_beg), int(ref_end)


def edlib_get_cigar(
    query: str,
    ref: str,
) -> Optional[str]:
    """
    HW alignment returning only CIGAR.
    """
    require_edlib()

    res = edlib.align(query, ref, mode="HW", task="path")
    cigar = res.get("cigar")
    return str(cigar) if cigar else None
