# botas/core/pairing.py
"""
Paired-end pairing logic (single source of truth).

All PE geometry, insert-size computation, and proper-pair
validation MUST live here.
"""

from __future__ import annotations
from typing import Tuple
from botas.core.circular import circular_insert_fr


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
    max_insert: int,
    L: int,
) -> Tuple[bool, int, str]:
    """
    Proper-pair check for circular references.

    ``expected_insert`` belongs to mate rescue: it predicts where an absent
    mate should be searched for.  Proper-pair validation must instead honor
    the user-facing ``max_insert`` limit, just as the linear path does.
    """
    ok, ins, orient = circular_insert_fr(
        pos1, strand1, pos2, strand2, read_len, L
    )

    if not ok:
        return False, 0, "invalid"

    if not 0 < ins <= max_insert:
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
) -> Tuple[bool, int, str]:
    """
    Unified PE proper-pair check (linear or circular).

    ``expected_insert`` remains in this interface for compatibility with the
    alignment pipeline, but validation intentionally uses ``max_insert``.
    """
    if circular:
        return is_proper_pair_circular(
            pos1,
            strand1,
            pos2,
            strand2,
            read_len,
            max_insert,
            ref_len,
        )

    return is_proper_pair_linear(
        pos1,
        strand1,
        pos2,
        strand2,
        read_len,
        max_insert,
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
