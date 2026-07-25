#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import logging

from botas.core.slice import slice_reference
from botas.core.align_core import Hit, align_read
from botas.core.utils import revcomp
from botas.core.circular import (
    circ_dist,
    expected_mate_start,
    check_pair_circular,
)
from botas.core.edlib_utils import edlib_align_hw_locations, edlib_get_cigar
from botas.core.pairing import is_proper_pair_unified

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PairHit:
    hit1: Optional[Hit]
    hit2: Optional[Hit]
    proper_pair: bool
    insert_size: Optional[int]


def _safe_pair_check(
    *,
    hit1: Hit,
    hit2: Hit,
    read_len: int,
    ref_len: int,
    circular: bool,
    max_insert: int,
    expected_insert: int,
):
    out = is_proper_pair_unified(
        pos1=hit1.pos0,
        strand1=hit1.strand,
        pos2=hit2.pos0,
        strand2=hit2.strand,
        read_len=read_len,
        ref_len=ref_len,
        circular=circular,
        max_insert=max_insert,
        expected_insert=expected_insert,
    )

    if out is None:
        return False, None, None

    return out


def rescue_mate(
    *,
    mate_seq: str,
    mate_strand: str,
    anchor_hit: Hit,
    ref_seq: str,
    index,
    circular: bool,
    expected_insert: int,
    rescue_pad: int = 200,
    use_fast: bool = True,
    original_len: int | None = None,
    circular_overhang: int = 0,
) -> Optional[Hit]:
    L_ref = len(ref_seq)
    L = int(original_len or L_ref)
    read_len = len(mate_seq)

    try:
        center0 = expected_mate_start(
            anchor_hit.pos0,
            anchor_hit.strand,
            mate_strand,
            expected_insert,
            read_len,
            L,
        )
    except ValueError:
        return None

    if use_fast:
        query = mate_seq if mate_strand == "+" else revcomp(mate_seq)

        best_p0 = None
        best_d = rescue_pad + 1

        for mm in index.iter_read_minimizers(query):
            p0 = mm[0]
            p0_orig = (p0 - circular_overhang) % L if circular else p0
            d = circ_dist(p0_orig, center0, L) if circular else abs(p0 - center0)

            if d <= rescue_pad and d < best_d:
                best_d = d
                best_p0 = p0

                if best_d == 0:
                    break

        if best_p0 is not None:
            center0 = best_p0

    query = mate_seq if mate_strand == "+" else revcomp(mate_seq)

    ref_slice, slice_start = slice_reference(
        ref_seq=ref_seq,
        center0=center0,
        read_len=len(query),
        pad=rescue_pad,
        circular=circular,
    )

    if not ref_slice:
        return None

    out = edlib_align_hw_locations(query, ref_slice)
    if out is None:
        return None

    edits, rb, _ = out

    max_edits = int(0.15 * len(query))
    if edits > max_edits:
        return None

    cigar = edlib_get_cigar(query, ref_slice)
    if cigar is None:
        return None

    abs_start = slice_start + rb

    if circular:
        pos0 = (abs_start - circular_overhang) % L
    else:
        pos0 = max(0, min(L - 1, abs_start))

    return Hit(
        rname=anchor_hit.rname,
        pos0=int(pos0),
        strand=mate_strand,
        cigar=cigar,
        ascore=-edits,
        mapq=min(60, 30 + max(0, 30 - edits)),
    )


def align_pair_simple(
    *,
    r1_seq: str,
    r2_seq: str,
    rname: str,
    ref_seq: str,
    index,
    circular: bool,
    k: int = 15,
    step: int = 5,
    pad: int = 250,
    max_windows: int = 50,
    min_seed_hits: int = 2,
    max_insert: int = 1200,
    expected_insert: int = 300,
    rescue_pad: int = 200,
    do_rescue: bool = True,
    sensitive: bool = False,
    original_len: int | None = None,
    circular_overhang: int = 0,
) -> PairHit:

    if sensitive:
        step = 1
        pad = max(pad, 500)
        max_windows = max(max_windows, 100)
        min_seed_hits = 1
        max_insert = max(max_insert, 5000)
        rescue_pad = max(rescue_pad, 400)

    L_ref = len(ref_seq)
    L = int(original_len or L_ref) if circular else L_ref
    read_len = len(r1_seq)

    # Existing BOTAS behavior: align reads linearly even when reference is unrolled circular.
    se_circular = False

    hit1 = align_read(
        read_seq=r1_seq,
        rname=rname,
        ref_seq=ref_seq,
        index=index,
        circular=se_circular,
        k=k,
        step=step,
        pad=pad,
        max_windows=max_windows,
        min_seed_hits=min_seed_hits,
    )

    hit2 = None

    hit2_rev = align_read(
        read_seq=revcomp(r2_seq),
        rname=rname,
        ref_seq=ref_seq,
        index=index,
        circular=se_circular,
        k=k,
        step=step,
        pad=pad,
        max_windows=max_windows,
        min_seed_hits=min_seed_hits,
    )

    if hit2_rev:
        hit2 = Hit(
            rname=hit2_rev.rname,
            pos0=hit2_rev.pos0,
            strand="-",
            cigar=hit2_rev.cigar,
            ascore=hit2_rev.ascore,
            mapq=hit2_rev.mapq,
        )

    def _norm_hit(hit):
        if hit is None or not circular:
            return hit
        return Hit(
            rname=hit.rname,
            pos0=(int(hit.pos0) - int(circular_overhang)) % L,
            strand=hit.strand,
            cigar=hit.cigar,
            ascore=hit.ascore,
            mapq=hit.mapq,
        )

    hit1 = _norm_hit(hit1)
    hit2 = _norm_hit(hit2)

    # Primary pair check
    if hit1 and hit2:
        proper, ins, _ = _safe_pair_check(
            hit1=hit1,
            hit2=hit2,
            read_len=read_len,
            ref_len=L,
            circular=circular,
            max_insert=max_insert,
            expected_insert=expected_insert,
        )

        if proper:
            return PairHit(
                hit1=hit1,
                hit2=hit2,
                proper_pair=True,
                insert_size=ins,
            )

    # Forced circular mate rescue
    if circular and hit1 and hit2:
        if hit1.ascore >= hit2.ascore:
            anchor, weak = hit1, hit2
            weak_is_r2 = True
        else:
            anchor, weak = hit2, hit1
            weak_is_r2 = False

        ANCHOR_SCORE_MIN = -5
        WEAK_OVERRIDE_DIST = 80

        orig_weak = weak
        invalidate = False

        if anchor.ascore >= ANCHOR_SCORE_MIN:
            d = circ_dist(anchor.pos0, weak.pos0, L)
            if d > WEAK_OVERRIDE_DIST:
                invalidate = True
                if weak_is_r2:
                    hit2 = None
                else:
                    hit1 = None

        rescued = rescue_mate(
            mate_seq=r2_seq if weak_is_r2 else r1_seq,
            mate_strand=orig_weak.strand,
            anchor_hit=anchor,
            ref_seq=ref_seq,
            index=index,
            circular=True,
            expected_insert=expected_insert,
            rescue_pad=max(rescue_pad, 500),
            use_fast=False,
            original_len=L,
            circular_overhang=circular_overhang,
        )

        if rescued:
            if weak_is_r2:
                hit2 = rescued
            else:
                hit1 = rescued
        elif invalidate:
            if weak_is_r2:
                hit2 = orig_weak
            else:
                hit1 = orig_weak

    # General mate rescue
    if do_rescue:
        anchor = hit1 or hit2

        if anchor:
            is_r1_anchor = anchor is hit1
            mate_seq = r2_seq if is_r1_anchor else r1_seq
            mate_strand = "-" if is_r1_anchor else "+"

            rescued = rescue_mate(
                mate_seq=mate_seq,
                mate_strand=mate_strand,
                anchor_hit=anchor,
                ref_seq=ref_seq,
                index=index,
                circular=circular,
                expected_insert=expected_insert,
                rescue_pad=rescue_pad,
                use_fast=True,
                original_len=L,
                circular_overhang=circular_overhang,
            )

            if rescued:
                if is_r1_anchor:
                    hit2 = rescued
                else:
                    hit1 = rescued

    proper = False
    ins = None

    if circular:
        if hit1 and hit1.pos0 >= L:
            hit1 = Hit(
                rname=hit1.rname,
                pos0=hit1.pos0 - L,
                strand=hit1.strand,
                cigar=hit1.cigar,
                ascore=hit1.ascore,
                mapq=hit1.mapq,
            )

        if hit2 and hit2.pos0 >= L:
            hit2 = Hit(
                rname=hit2.rname,
                pos0=hit2.pos0 - L,
                strand=hit2.strand,
                cigar=hit2.cigar,
                ascore=hit2.ascore,
                mapq=hit2.mapq,
            )

    if hit1 and hit2:
        if circular:
            proper, ins, _ = check_pair_circular(
                hit1.pos0,
                hit1.strand,
                hit2.pos0,
                hit2.strand,
                read_len,
                expected_insert,
                L,
                tol_ins=50,
            )
        else:
            proper, ins, _ = _safe_pair_check(
                hit1=hit1,
                hit2=hit2,
                read_len=read_len,
                ref_len=L,
                circular=False,
                max_insert=max_insert,
                expected_insert=expected_insert,
            )

    return PairHit(
        hit1=hit1,
        hit2=hit2,
        proper_pair=proper,
        insert_size=ins,
    )


def align_pair(
    *,
    r1_seq: str,
    r2_seq: str,
    rname: str,
    ref_seq: str,
    index,
    circular: bool,
    k: int = 15,
    step: int = 5,
    pad: int = 250,
    max_windows: int = 50,
    min_seed_hits: int = 2,
    max_insert: int = 1200,
    expected_insert: int = 300,
    rescue_pad: int = 200,
    do_rescue: bool = True,
    sensitive: bool = False,
    debug: bool = False,
    pair_id: int = 0,
    rescue_strategy: str = "fast_first",
    min_rescue_seeds: int = 2,
    max_rescue_edits: int = 15,
    fast_rescue_only: bool = False,
    **kwargs,
) -> PairHit:

    if debug:
        return align_pair_debug(
            r1_seq=r1_seq,
            r2_seq=r2_seq,
            rname=rname,
            ref_seq=ref_seq,
            index=index,
            circular=circular,
            k=k,
            step=step,
            pad=pad,
            max_windows=max_windows,
            min_seed_hits=min_seed_hits,
            max_insert=max_insert,
            expected_insert=expected_insert,
            rescue_pad=rescue_pad,
            do_rescue=do_rescue,
            sensitive=sensitive,
            pair_id=pair_id,
            original_len=kwargs.get("original_len"),
            circular_overhang=kwargs.get("circular_overhang", 0),
        )

    return align_pair_simple(
        r1_seq=r1_seq,
        r2_seq=r2_seq,
        rname=rname,
        ref_seq=ref_seq,
        index=index,
        circular=circular,
        k=k,
        step=step,
        pad=pad,
        max_windows=max_windows,
        min_seed_hits=min_seed_hits,
        max_insert=max_insert,
        expected_insert=expected_insert,
        rescue_pad=rescue_pad,
        do_rescue=do_rescue,
        sensitive=sensitive,
        original_len=kwargs.get("original_len"),
        circular_overhang=kwargs.get("circular_overhang", 0),
    )


def align_pair_sensitive(
    *,
    r1_seq: str,
    r2_seq: str,
    rname: str,
    ref_seq: str,
    index,
    circular: bool,
    **kwargs,
) -> PairHit:
    kwargs["sensitive"] = True

    return align_pair(
        r1_seq=r1_seq,
        r2_seq=r2_seq,
        rname=rname,
        ref_seq=ref_seq,
        index=index,
        circular=circular,
        **kwargs,
    )


def align_pair_debug(
    *,
    r1_seq: str,
    r2_seq: str,
    rname: str,
    ref_seq: str,
    index,
    circular: bool,
    pair_id: int = 0,
    **kwargs,
) -> PairHit:

    old_level = logger.level
    logger.setLevel(logging.DEBUG)

    logger.debug("\n=== DEBUG Pair %s ===", pair_id)
    logger.debug("R1 length: %d, R2 length: %d", len(r1_seq), len(r2_seq))

    try:
        result = align_pair_simple(
            r1_seq=r1_seq,
            r2_seq=r2_seq,
            rname=rname,
            ref_seq=ref_seq,
            index=index,
            circular=circular,
            **kwargs,
        )

        logger.debug("Result: proper=%s, insert=%s", result.proper_pair, result.insert_size)

        if result.hit1:
            logger.debug(
                "R1: pos=%s, strand=%s, mapq=%s",
                result.hit1.pos0,
                result.hit1.strand,
                result.hit1.mapq,
            )
        else:
            logger.debug("R1: unmapped")

        if result.hit2:
            logger.debug(
                "R2: pos=%s, strand=%s, mapq=%s",
                result.hit2.pos0,
                result.hit2.strand,
                result.hit2.mapq,
            )
        else:
            logger.debug("R2: unmapped")

        return result

    finally:
        logger.setLevel(old_level)