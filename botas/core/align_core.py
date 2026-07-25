"""
Core alignment module for Botas - seed-and-extend aligner using Edlib.

IMPORTANT DESIGN RULE:
- Single-end alignment is performed in LINEAR reference coordinates only.
- Circularity must be enforced at the FRAGMENT/PAIR level (PE) using botas.core.circular.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple
import logging
import os
import sys

from botas.core.ref_index import KmerIndex, Window, merge_windows
from botas.core.scoring import score_to_mapq
from botas.core.utils import revcomp
from botas.core.cigar import ref_leading_consumption
from botas.core.edlib_utils import edlib_align_hw_path
from botas.core.slice import slice_reference
from botas.core.circular import circ_dist

logger = logging.getLogger(__name__)

_PROFILE_ALIGN = os.environ.get("BOTAS_PROFILE_ALIGN", "0") == "1"
_PROF_READS = 0
_PROF_WINDOWS = 0
_PROF_EDLIB = 0


@dataclass(frozen=True)
class Hit:
    rname: str
    pos0: int
    strand: str
    cigar: str
    ascore: int
    mapq: int


@dataclass(frozen=True)
class _Candidate:
    rname: str
    slice_start: int
    strand: str
    score: int
    cigar: str
    ref_beg: int
    ref_end: int
    win_hits: int
    win_start0: int


def _profile_report() -> None:
    if not _PROFILE_ALIGN:
        return

    avg_w = _PROF_WINDOWS / max(1, _PROF_READS)
    avg_e = _PROF_EDLIB / max(1, _PROF_READS)

    print(
        f"[align-profile] reads={_PROF_READS} "
        f"windows={_PROF_WINDOWS} "
        f"edlib={_PROF_EDLIB} "
        f"avg_windows={avg_w:.2f} "
        f"avg_edlib={avg_e:.2f}",
        file=sys.stderr,
        flush=True,
    )


def _align_one_window(win, *, query_seq, rname, ref_seq, L, circular, strand):
    global _PROF_EDLIB

    if circular:
        length = int(win.end0 - win.start0)
        if win.start0 < 0 or win.end0 > L:
            length += 250
        if length <= 0:
            return None

        ref_slice, slice_start = slice_reference(
            ref_seq=ref_seq,
            start0=int(win.start0),
            length=length,
            circular=True,
        )
        if not ref_slice:
            return None
    else:
        slice_start = max(0, int(win.start0))
        slice_end = min(L, int(win.end0))
        if slice_end <= slice_start:
            return None

        ref_slice = ref_seq[slice_start:slice_end]
        if not ref_slice:
            return None

    if _PROFILE_ALIGN:
        _PROF_EDLIB += 1

    out = edlib_align_hw_path(query_seq, ref_slice)
    if out is None:
        return None

    score, cigar, ref_beg, ref_end = out

    cand = _Candidate(
        rname=rname,
        slice_start=slice_start,
        strand=strand,
        score=int(score),
        cigar=str(cigar),
        ref_beg=int(ref_beg),
        ref_end=int(ref_end),
        win_hits=int(getattr(win, "hits", 0)),
        win_start0=int(win.start0),
    )

    lead = ref_leading_consumption(cand.cigar)
    implied_abs_start = cand.slice_start + cand.ref_beg + lead
    win_start = cand.win_start0

    if circular:
        implied0 = implied_abs_start % L
        win0 = win_start % L
        consistency_penalty = circ_dist(implied0, win0, L)
    else:
        win_start = max(0, min(L, win_start))
        consistency_penalty = abs(implied_abs_start - win_start)

    cur_tie = (cand.score, -cand.win_hits, int(consistency_penalty))
    return cand, cur_tie


def align_read(
    *,
    read_seq: str,
    rname: str,
    ref_seq: str,
    index: KmerIndex,
    circular: bool,
    k: int = 15,
    step: int = 5,
    pad: int = 250,
    max_windows: int = 50,
    min_seed_hits: int = 2,
    threads: int = 1,
    strands: tuple[str, ...] = ("+", "-"),
) -> Optional[Hit]:

    global _PROF_READS, _PROF_WINDOWS

    if not read_seq:
        return None

    L = len(ref_seq)
    if L == 0:
        return None

    best: Optional[_Candidate] = None
    best_tie: Optional[Tuple[int, int, int]] = None

    strand_queries = []

    if "+" in strands:
        strand_queries.append(("+", read_seq))

    if "-" in strands:
        strand_queries.append(("-", revcomp(read_seq)))

    for strand, query_seq in strand_queries:
        wins = index.windows_for_read(
            query_seq,
            pad=pad,
            k=k,
            step=step,
            min_hits=min_seed_hits,
        )

        if _PROFILE_ALIGN:
            _PROF_READS += 1
            _PROF_WINDOWS += len(wins)

            if _PROF_READS % 10000 == 0:
                _profile_report()

        if not wins:
            continue

        wins = merge_windows(wins)

        if circular:
            extra = []
            for w in wins:
                if w.start0 < pad:
                    extra.append(Window(w.start0 + L, w.end0 + L))
                if w.end0 > (L - pad):
                    extra.append(Window(w.start0 - L, w.end0 - L))
            if extra:
                wins = merge_windows(wins + extra)

        if len(wins) > max_windows:
            wins = wins[:max_windows]

        for win in wins:
            res = _align_one_window(
                win,
                query_seq=query_seq,
                rname=rname,
                ref_seq=ref_seq,
                L=L,
                circular=circular,
                strand=strand,
            )

            if res is None:
                continue

            cand, cur_tie = res

            if best is None or cand.score > best.score:
                best = cand
                best_tie = cur_tie
            elif cand.score == best.score:
                if best_tie is None or cur_tie < best_tie:
                    best = cand
                    best_tie = cur_tie

    if best is None:
        return None

    lead = ref_leading_consumption(best.cigar)
    pos0 = best.slice_start + best.ref_beg + lead

    if circular:
        pos0 %= L
    else:
        if pos0 < 0 or pos0 >= L:
            return None

    mq = score_to_mapq(
        best.score,
        read_len=len(read_seq),
        num_windows=1,
    )

    return Hit(
        rname=rname,
        pos0=int(pos0),
        strand=best.strand,
        cigar=best.cigar,
        ascore=int(best.score),
        mapq=int(mq),
    )


def align_read_batch(
    reads: List[Tuple[str, str]],
    rname: str,
    ref_seq: str,
    index: KmerIndex,
    circular: bool,
    **kwargs,
) -> List[Tuple[str, Optional[Hit]]]:
    results: List[Tuple[str, Optional[Hit]]] = []

    for read_id, read_seq in reads:
        hit = align_read(
            read_seq=read_seq,
            rname=rname,
            ref_seq=ref_seq,
            index=index,
            circular=circular,
            **kwargs,
        )
        results.append((read_id, hit))

    return results