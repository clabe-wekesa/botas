# botas/core/pe_pool.py
from __future__ import annotations
from multiprocessing import Pool
from dataclasses import dataclass
import logging
from botas.io.reference_set import load_reference_set
from botas.core.ref_index import KmerIndex
from botas.core.align_pe import align_pair
from botas.core.align_core import Hit
from botas.core.parallel_context import PEContext
from botas.core.utils import chunked

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _ContigCand:
    contig: object
    h1: int
    h2: int
    score: int


def _contig_candidates(refset, r1: str, r2: str, min_seed_hits: int, top_n: int = 3):
    """
    Rank contigs by JOINT seed evidence from both reads.
    Only contigs with evidence from both reads are prioritized.
    """
    cands = []
    for c in refset.contigs():
        try:
            h1 = c.index.num_seed_hits(r1)
            h2 = c.index.num_seed_hits(r2)
        except AttributeError:
            return list(refset.contigs())

        if (h1 + h2) < min_seed_hits:
            continue

        score = h1 * h2
        cands.append(_ContigCand(c, h1, h2, score))

    if not cands:
        return list(refset.contigs())

    cands.sort(key=lambda x: (x.score, x.h1 + x.h2), reverse=True)
    return [x.contig for x in cands[:top_n]]


def _has_seed_evidence(index, r1: str, r2: str, min_seed_hits: int) -> bool:
    try:
        n1 = index.num_seed_hits(r1)
        n2 = index.num_seed_hits(r2)
    except AttributeError:
        return True
    return (n1 + n2) >= min_seed_hits


@dataclass(frozen=True)
class _PEWorkerConfig:
    ref_fa: str
    circular_all: bool
    circular_overhang_percent: int
    k: int
    step: int
    pad: int
    max_windows: int
    min_seed_hits: int
    max_insert: int
    expected_insert: int
    rescue_pad: int
    do_rescue: bool
    sensitive: bool
    debug_pairs: int


# ============================================================
# Pre-built contig data (passed from main process to workers)
# ============================================================

@dataclass
class _ContigData:
    """Serializable contig snapshot: seq + pre-built index."""
    name: str
    seq: str
    circular: bool
    length: int
    index: KmerIndex
    orig_len: int | None = None
    circular_overhang: int = 0


# ============================================================
# Worker context (GLOBAL per process)
# ============================================================

_CTX = None  # (cfg, list[_ContigData])


def _init_worker(cfg: _PEWorkerConfig, contigs: list[_ContigData]):
    """
    Per-process initialization.

    FIX: Previously this function called load_reference_set() and rebuilt
    a KmerIndex for every contig inside every worker process — meaning N
    worker processes each did the full index build independently.

    Now the main process builds all indexes once (in align_paired_pool)
    and passes the pre-built _ContigData objects to workers via the Pool
    initializer. Workers receive fully-built indexes at startup and do no
    index construction at all.
    """
    global _CTX
    _CTX = (cfg, contigs)


def _normalize_circular_hit(hit, contig):
    """Map padded-reference positions back to original circular coordinates.

    This is intentionally cheap: one subtraction and one modulo, only for
    final reported hits, not during candidate generation or DP.
    """
    if hit is None or not getattr(contig, "circular", False):
        return hit

    orig_len = getattr(contig, "orig_len", None) or getattr(contig, "length", None) or len(contig.seq)
    overhang = getattr(contig, "circular_overhang", 0) or 0

    pos0 = (int(hit.pos0) - int(overhang)) % int(orig_len)

    return Hit(
        rname=hit.rname,
        pos0=pos0,
        strand=hit.strand,
        cigar=hit.cigar,
        ascore=hit.ascore,
        mapq=hit.mapq,
    )


def _pick_better_pair(best, cand):
    if best is None:
        return cand

    def _rank(ph):
        mapped = int(ph.hit1 is not None) + int(ph.hit2 is not None)
        proper = 1 if ph.proper_pair else 0
        as_sum = (ph.hit1.ascore if ph.hit1 else -(10**9)) + (ph.hit2.ascore if ph.hit2 else -(10**9))
        mq_sum = (ph.hit1.mapq if ph.hit1 else 0) + (ph.hit2.mapq if ph.hit2 else 0)
        return (proper, mapped, as_sum, mq_sum)

    return cand if _rank(cand) > _rank(best) else best


def _pe_worker(chunk):
    cfg, contigs = _CTX
    out = []

    for idx, qname, r1s, r1q, r2s, r2q in chunk:
        best_ph = None

        # Fast path for single-contig references.
        # For E. coli there is only one contig, so contig pre-filtering
        # and ranking are unnecessary and only add extra seed-counting cost.
        if len(contigs) == 1:
            cand_contigs = contigs
        else:
            scored = []

            for c in contigs:
                try:
                    h1 = c.index.num_seed_hits(r1s)
                    h2 = c.index.num_seed_hits(r2s)
                except AttributeError:
                    scored.append((c, 0, 0))
                    continue

                if (h1 + h2) < cfg.min_seed_hits:
                    continue

                scored.append((c, h1 * h2, h1 + h2))

            if scored:
                scored.sort(key=lambda x: (x[1], x[2]), reverse=True)
                cand_contigs = [x[0] for x in scored[:3]]
            else:
                cand_contigs = contigs

        for c in cand_contigs:
            ph = align_pair(
                r1_seq=r1s,
                r2_seq=r2s,
                rname=c.name,
                ref_seq=c.seq,
                index=c.index,
                circular=c.circular,
                k=getattr(c.index, "k", cfg.k),
                step=cfg.step,
                pad=cfg.pad,
                max_windows=cfg.max_windows,
                min_seed_hits=cfg.min_seed_hits,
                max_insert=cfg.max_insert,
                expected_insert=cfg.expected_insert,
                rescue_pad=cfg.rescue_pad,
                do_rescue=cfg.do_rescue,
                sensitive=cfg.sensitive,
                debug=(cfg.debug_pairs > 0 and idx < cfg.debug_pairs),
                pair_id=idx,
                original_len=getattr(c, "orig_len", None),
                circular_overhang=getattr(c, "circular_overhang", 0),
            )

            best_ph = _pick_better_pair(best_ph, ph)

            if best_ph and best_ph.proper_pair and best_ph.hit1 and best_ph.hit2:
                break

        out.append((idx, qname, r1s, r1q, r2s, r2q, best_ph))

    return out


def align_paired_pool(*, pairs_iter, ctx: PEContext, threads: int, chunk_size: int = 5000):
    """
    pairs_iter yields: (idx, r1, r2)

    FIX: KmerIndex is now built ONCE here in the main process, then passed
    to all worker processes via the Pool initializer. Previously each of
    the N worker processes independently read the FASTA from disk and
    rebuilt the full index, which was wasteful for multi-contig references.
    """
    cfg = _PEWorkerConfig(
        ref_fa=ctx.ref_fa,
        circular_all=ctx.circular_all,
        circular_overhang_percent=ctx.circular_overhang_percent,
        k=ctx.k,
        step=ctx.step,
        pad=ctx.pad,
        max_windows=ctx.max_windows,
        min_seed_hits=ctx.min_seed_hits,
        max_insert=ctx.max_insert,
        expected_insert=ctx.expected_insert,
        rescue_pad=ctx.rescue_pad,
        do_rescue=ctx.do_rescue,
        sensitive=ctx.sensitive,
        debug_pairs=ctx.debug_pairs,
    )

    # Build indexes once in the main process
    if ctx.prebuilt_contigs is not None:
        logger.info("Using prebuilt BOTAS index contigs.")
        contigs = list(ctx.prebuilt_contigs)

        if ctx.circular_all:
            for c in contigs:
                c.circular = True
    else:
        logger.info("Building reference indexes (main process)...")
        refset = load_reference_set(
            ctx.ref_fa,
            circular=ctx.circular_all,
            circular_overhang_percent=ctx.circular_overhang_percent,
        )

        contigs: list[_ContigData] = []
        for c in refset.contigs():
            idx = KmerIndex(c.seq, k=ctx.k, circular=c.circular, label=c.name)
            contigs.append(_ContigData(
                name=c.name,
                seq=c.seq,
                circular=c.circular,
                length=c.length,
                index=idx,
                orig_len=getattr(c, "orig_len", c.length),
                circular_overhang=getattr(c, "circular_overhang", 0),
            ))

    logger.info("Indexes ready for %d contig(s). Spawning %d worker(s).", len(contigs), threads)

    pool = Pool(
        processes=threads,
        initializer=_init_worker,
        initargs=(cfg, contigs),
    )

    try:
        for results in pool.imap_unordered(
            _pe_worker,
            (
                [
                    (
                        idx,
                        r1.name,
                        r1.seq,
                        r1.qual,
                        r2.seq,
                        r2.qual,
                    )
                    for idx, r1, r2 in chunk
                ]
                for chunk in chunked(pairs_iter, chunk_size)
            ),
            chunksize=1,
        ):
            for rec in results:
                yield rec
    finally:
        pool.close()
        pool.join()
