# botas/core/circular.py
from __future__ import annotations
from typing import Tuple
import math

def circ_mean(pos, L):
    x = 0.0
    y = 0.0
    for p in pos:
        a = 2.0 * math.pi * (p / L)
        x += math.cos(a)
        y += math.sin(a)
    return int((math.atan2(y, x) / (2.0 * math.pi) * L) % L)

def norm_pos(pos0: int, L: int) -> int:
    """Normalize 0-based position into [0, L-1]."""
    if L <= 0:
        raise ValueError("L must be > 0")
    return pos0 % L


def forward_delta(a: int, b: int, L: int) -> int:
    """Forward distance from a to b along + direction."""
    return (norm_pos(b, L) - norm_pos(a, L)) % L


def circular_dist(a: int, b: int, L: int) -> int:
    """Shortest distance between two positions on a circle."""
    d = forward_delta(a, b, L)
    return min(d, L - d)


# -----------------------------
# Core circular geometry
# -----------------------------

def circ_dist(a: int, b: int, L: int) -> int:
    return circular_dist(a, b, L)

def circ_forward_delta(a: int, b: int, L: int) -> int:
    return forward_delta(a, b, L)

def circ_backward_delta(a: int, b: int, L: int) -> int:
    """Backward distance from a to b along - direction (mod L)."""
    return (a - b) % L

def circular_insert_fr(
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    L: int,
) -> Tuple[bool, int, str]:
    """
    Canonical circular insert computation.
    Returns (is_proper, insert_size, orientation).
    """
    if strand1 == "+" and strand2 == "-":
        return True, forward_delta(pos1, pos2, L) + read_len, "FR"
    if strand1 == "-" and strand2 == "+":
        return True, forward_delta(pos2, pos1, L) + read_len, "RF"
    return False, 0, "invalid"


# -----------------------------
# Insert geometry (PE)
# -----------------------------

def infer_insert_FR(pos1: int, pos2: int, read_len: int, L: int) -> int:
    return forward_delta(pos1, pos2, L) + read_len


def infer_insert_RF(pos1: int, pos2: int, read_len: int, L: int) -> int:
    return forward_delta(pos2, pos1, L) + read_len

# -----------------------------
# Expected mate start
# -----------------------------

def expected_mate_start(
    anchor_pos: int,
    anchor_strand: str,
    mate_strand: str,
    insert: int,
    read_len: int,
    L: int,
) -> int:
    """
    Expected mate start position given anchor start.
    Matches truth definition:
      mate_start = anchor + (insert - read_len)
    """
    shift = insert - read_len

    if anchor_strand == "+" and mate_strand == "-":
        return norm_pos(anchor_pos + shift, L)
    if anchor_strand == "-" and mate_strand == "+":
        return (anchor_pos - shift) % L

    raise ValueError("Invalid strand combination")

# -----------------------------
# Pair validation
# -----------------------------

def check_pair_circular(
    pos1: int,
    strand1: str,
    pos2: int,
    strand2: str,
    read_len: int,
    insert: int,
    L: int,
    tol_ins: int = 20,
) -> Tuple[bool, int, str]:
    """
    Validate and score a circular PE pair.

    Returns:
        (is_proper, observed_insert, orientation)
    """
    # Delegate geometry to the canonical function
    ok, ins_obs, orientation = circular_insert_fr(
        pos1,
        strand1,
        pos2,
        strand2,
        read_len,
        L,
    )

    if not ok:
        return False, 0, "invalid"

    # Insert-size validation
    if abs(ins_obs - insert) > tol_ins:
        return False, ins_obs, orientation

    return True, ins_obs, orientation

