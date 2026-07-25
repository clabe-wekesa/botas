# botas/core/se_pool.py
from __future__ import annotations
from multiprocessing import Pool
from dataclasses import dataclass
import logging
from botas.io.reference_set import load_reference_set
from botas.core.ref_index import KmerIndex
from botas.core.align_core import align_read
from botas.core.utils import chunked

logger = logging.getLogger(__name__)


# ============================================================
# Worker config + context (mirrors PE)
# ============================================================

@dataclass(frozen=True)
class _SEWorkerConfig:
    k: int
    step: int
    pad: int
    max_windows: int
    min_seed_hits: int


@dataclass
class _ContigData:
    """Serializable contig snapshot: seq + pre-built index."""
    name: str
    seq: str
    circular: bool
    length: int
    index: KmerIndex


_CTX = None  # (cfg, list[_ContigData])


def _init_worker(cfg: _SEWorkerConfig, contigs: list[_ContigData]):
    """
    Per-process initialization.

    FIX: Previously loaded the reference FASTA and rebuilt KmerIndex for
    every contig inside each worker. Now receives pre-built contig data
    (including indexes) built once in the main process via align_single_pool.
    """
    global _CTX
    _CTX = (cfg, contigs)

def _seed_hit_count(seq: str, index: KmerIndex, k: int, step: int) -> int:
    """
    Cheaply estimate whether a read belongs to a contig.

    This only counts seed hits. It does not run Edlib.
    """
    total = 0

    if len(seq) < k:
        return 0

    table = getattr(index, "index", None)
    if table is None:
        table = getattr(index, "idx", None)
    if table is None:
        table = getattr(index, "_index", None)

    if table is None:
        raise AttributeError(
            "KmerIndex has no accessible seed table: expected index, idx, or _index"
        )

    for i in range(0, len(seq) - k + 1, step):
        seed = seq[i:i + k]
        hits = table.get(seed)
        if hits:
            total += len(hits)

    return total


def _rank_contigs_for_read(seq: str, contigs: list[_ContigData], k: int, step: int, top_n: int = 2):
    """
    Rank contigs by seed support before expensive alignment.

    Returns only top_n contigs with seed evidence.
    """
    scored = []

    for c in contigs:
        score = _seed_hit_count(seq, c.index, k, step)
        if score > 0:
            scored.append((score, c))

    if not scored:
        return contigs

    scored.sort(key=lambda x: x[0], reverse=True)
    return [scored[0][1]]    


# ============================================================
# Worker
# ============================================================

def _se_worker(chunk):
    """
    chunk: list of (idx, qname, seq, qual)
    """
    cfg, contigs = _CTX
    out = []

    for idx, qname, seq, qual in chunk:
        best_hit = None
        second_hit = None

        candidate_contigs = _rank_contigs_for_read(
            seq,
            contigs,
            cfg.k,
            cfg.step,
            top_n=2,
        )

        for c in candidate_contigs:
            hit = align_read(
                read_seq=seq,
                rname=c.name,
                ref_seq=c.seq,
                index=c.index,
                circular=False,
                k=cfg.k,
                step=cfg.step,
                pad=cfg.pad,
                max_windows=cfg.max_windows,
                min_seed_hits=cfg.min_seed_hits,
            )

            if hit is None:
                continue

            if best_hit is None or hit.ascore > best_hit.ascore:
                second_hit = best_hit
                best_hit = hit
            elif second_hit is None or hit.ascore > second_hit.ascore:
                second_hit = hit

        out.append((idx, qname, seq, qual, best_hit, second_hit))

    return out


# ============================================================
# Public API
# ============================================================

def align_single_pool(
    *,
    reads_iter,
    ref_fa: str,
    threads: int,
    chunk_size: int,
    k: int,
    step: int,
    pad: int,
    max_windows: int,
    min_seed_hits: int,
    prebuilt_contigs=None,
):
    """
    Parallel single-end alignment.

    reads_iter yields: (idx, read)

    FIX: KmerIndex is now built ONCE here in the main process and sent to
    all worker processes via the Pool initializer. Previously each worker
    independently read the FASTA from disk and rebuilt every contig index,
    wasting time proportional to (num_workers × genome_size).
    """
    cfg = _SEWorkerConfig(k=k, step=step, pad=pad, max_windows=max_windows, min_seed_hits=min_seed_hits)

    # Build indexes once in the main process
    contigs = []

    if prebuilt_contigs is not None:
        logger.info("Using prebuilt BOTAS index contigs for SE.")

        for c in prebuilt_contigs:
            contigs.append(
                _ContigData(
                    name=c.name,
                    seq=c.seq,
                    circular=getattr(c, "circular", False),
                    length=getattr(c, "length", len(c.seq)),
                    index=c.index,
                )
            )

    else:
        logger.info("Building SE reference indexes (main process)...")

        refset = load_reference_set(ref_fa)

        for c in refset.contigs():
            idx = KmerIndex(
                c.seq,
                k=k,
                circular=False,
                label=c.name,
            )

            contigs.append(
                _ContigData(
                    name=c.name,
                    seq=c.seq,
                    circular=False,
                    length=c.length,
                    index=idx,
                )
            )
        contigs.append(_ContigData(
            name=c.name,
            seq=c.seq,
            circular=False,
            length=c.length,
            index=idx,
        ))
    #logger.info("SE indexes built for %d contig(s). Spawning %d worker(s).", len(contigs), threads)

    pool = Pool(processes=threads, initializer=_init_worker, initargs=(cfg, contigs))
    try:
        for results in pool.imap_unordered(
            _se_worker,
            (
                [
                    (idx, r.name, r.seq, r.qual)
                    for idx, r in chunk
                ]
                for chunk in chunked(reads_iter, chunk_size)
            ),
            chunksize=1,
        ):
            for rec in results:
                yield rec
    finally:
        pool.close()
        pool.join()
