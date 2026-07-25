# botas/core/pairing.py
"""
Paired-end pairing logic (single source of truth).

All PE geometry, insert-size computation, and proper-pair
validation MUST live here.
"""

from __future__ import annotations
from typing import Tuple
from botas.core.circular import circular_insert_fr, circular_dist, circ_dist


def _best_shift(p_fixed, p_var, L, expected_insert):
    """Return p_var shifted by {0, +L, -L} that best matches expected insert."""
    best_p = p_var
    best_err = abs(circ_dist(p_fixed, p_var, L) - expected_insert)

    for shift in (L, -L):
        p_try = p_var + shift
        err = abs(abs(p_try - p_fixed) - expected_insert)
        if err < best_err:
            best_err = err
            best_p = p_try

    return best_p


# -------------------------------------------------------
# Linear pairing
# -------------------------------------------------------

def is_proper_pair_linear(
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    max_insert: int,
) -> Tuple[bool, int, str]:
    """
    Proper-pair check for linear references.

    Returns:
        (is_proper, insert_size, orientation)
    """
    # FR
    if strand1 == "+" and strand2 == "-" and pos1 <= pos2:
        ins = (pos2 - pos1) + read_len
        if 0 < ins <= max_insert:
            return True, ins, "FR"

    # RF
    if strand1 == "-" and strand2 == "+" and pos2 <= pos1:
        ins = (pos1 - pos2) + read_len
        if 0 < ins <= max_insert:
            return True, ins, "RF"

    return False, 0, "invalid"


# -------------------------------------------------------
# Circular pairing
# -------------------------------------------------------

def is_proper_pair_circular(
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    expected_insert: int,
    L: int,
    tol_ins: int,
) -> Tuple[bool, int, str]:
    """
    Proper-pair check for circular references.
    """
    ok, ins, orient = circular_insert_fr(
        pos1, strand1, pos2, strand2, read_len, L
    )

    if not ok:
        return False, 0, "invalid"

    if abs(ins - expected_insert) > tol_ins:
        return False, ins, orient

    return True, ins, orient


# -------------------------------------------------------
# Unified interface
# -------------------------------------------------------

def is_proper_pair_unified(
    *,
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    circular: bool,
    ref_len: int,
    max_insert: int,
    expected_insert: int,
    tol_ins: int = 50,
) -> Tuple[bool, int, str]:
    """
    Unified PE proper-pair check (linear or circular).
    """
    if circular:
        L = ref_len

        # ---- PE circular coordinate normalization ----
        # Try adjusting pos2 relative to pos1
        p2_adj = _best_shift(pos1, pos2, L, expected_insert)
        err2 = abs(abs(p2_adj - pos1) - expected_insert)

        # Try adjusting pos1 relative to pos2
        p1_adj = _best_shift(pos2, pos1, L, expected_insert)
        err1 = abs(abs(pos2 - p1_adj) - expected_insert)

        if err2 <= err1 and err2 <= tol_ins:
            pos2 = p2_adj
        elif err1 < err2 and err1 <= tol_ins:
            pos1 = p1_adj

        return is_proper_pair_circular(
            pos1,
            strand1,
            pos2,
            strand2,
            read_len,
            expected_insert,
            L,
            tol_ins,
        )


def compute_insert_size(
    *,
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    circular: bool,
    ref_len: int,
) -> int:
    """
    Compute insert size without enforcing thresholds.
    Returns 0 if orientation invalid.
    """
    if circular:
        ok, ins, _ = circular_insert_fr(
            pos1, strand1, pos2, strand2, read_len, ref_len
        )
        return ins if ok else 0

    # Linear
    if strand1 == "+" and strand2 == "-" and pos1 <= pos2:
        return (pos2 - pos1) + read_len
    if strand1 == "-" and strand2 == "+" and pos2 <= pos1:
        return (pos1 - pos2) + read_len

    return 0
