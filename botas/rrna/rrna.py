#!/usr/bin/env python3
"""
rrna.py — rRNA detection helpers for BOTAS.

The classifier uses exact k-mer membership against an rRNA reference database.
Database loading is fail-fast, ambiguous k-mers are ignored, and reads are
checked in both orientations without duplicating the (potentially very large)
rRNA k-mer database in memory.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from Bio import SeqIO

_DNA = frozenset("ACGT")
_RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of an uppercase DNA sequence."""
    return seq.translate(_RC_TABLE)[::-1]


def _validate_params(k: int, min_hits: int | None = None) -> None:
    if k <= 0:
        raise ValueError(f"k must be > 0; got {k}")
    if min_hits is not None and min_hits <= 0:
        raise ValueError(f"min_hits must be > 0; got {min_hits}")


def load_rrna_kmers(rrna_fa: str, k: int = 18) -> set[str]:
    """
    Build a set of valid A/C/G/T k-mers from an rRNA FASTA file.

    Parameters
    ----------
    rrna_fa
        rRNA reference FASTA path. Plain text and ``.gz`` are supported.
    k
        k-mer size.

    Returns
    -------
    set[str]
        Unique forward-orientation rRNA k-mers. Read sequences are tested in
        both orientations by :func:`is_rrna_like`, avoiding a doubled database
        in memory.

    Raises
    ------
    FileNotFoundError
        If the database path does not exist.
    ValueError
        If no usable k-mers are produced.
    RuntimeError
        If the database cannot be parsed/read.
    """
    _validate_params(k)

    path = Path(rrna_fa)
    if not path.is_file():
        raise FileNotFoundError(f"rRNA database not found: {rrna_fa}")

    logging.info("[rrna] loading rRNA database → %s", rrna_fa)

    rrna_kmers: set[str] = set()
    records = 0
    skipped_ambiguous = 0

    opener = gzip.open if str(path).lower().endswith(".gz") else open

    try:
        with opener(path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                records += 1
                seq = str(record.seq).upper().replace("U", "T")

                if len(seq) < k:
                    continue

                for i in range(len(seq) - k + 1):
                    kmer = seq[i : i + k]
                    if set(kmer) <= _DNA:
                        rrna_kmers.add(kmer)
                    else:
                        skipped_ambiguous += 1

    except Exception as exc:
        # Never silently continue with an empty or partial database: doing so
        # turns every read into a false negative while the pipeline appears to
        # have succeeded.
        raise RuntimeError(f"failed to load rRNA database {rrna_fa}: {exc}") from exc

    if records == 0:
        raise ValueError(f"rRNA database contains no FASTA records: {rrna_fa}")

    if not rrna_kmers:
        raise ValueError(
            f"rRNA database produced no valid A/C/G/T {k}-mers: {rrna_fa}"
        )

    logging.info(
        "[rrna] built k-mer set (k=%d) size=%s records=%s skipped_ambiguous=%s",
        k,
        f"{len(rrna_kmers):,}",
        f"{records:,}",
        f"{skipped_ambiguous:,}",
    )

    return rrna_kmers


def _has_min_hits(seq: str, rrna_kmers: set[str], k: int, min_hits: int) -> bool:
    """Return True once *min_hits* database k-mers have been observed."""
    if len(seq) < k:
        return False

    hits = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]

        # Reads containing N/other ambiguity codes should not generate matches
        # for that window. Avoid allocating set(kmer) in this hot loop.
        if any(base not in "ACGT" for base in kmer):
            continue

        if kmer in rrna_kmers:
            hits += 1
            if hits >= min_hits:
                return True

    return False


def is_rrna_like(
    seq: str,
    rrna_kmers: set[str],
    k: int = 18,
    min_hits: int = 10,
) -> bool:
    """
    Classify a read as rRNA-like using exact k-mer support.

    A read is classified as rRNA-like when at least ``min_hits`` of its sliding
    k-mers occur in the rRNA database. Both sequence orientations are tested so
    classification does not depend on library/read orientation.

    Notes
    -----
    ``min_hits=10`` is the recommended BOTAS default for 150-bp bacterial
    RNA-seq based on real-data calibration. It should still be benchmarked on
    additional read lengths/datasets before being treated as universal.
    """
    _validate_params(k, min_hits)

    if not rrna_kmers:
        raise ValueError("rrna_kmers is empty; refusing to classify reads")

    seq = seq.upper().replace("U", "T")

    if _has_min_hits(seq, rrna_kmers, k, min_hits):
        return True

    rc = _reverse_complement(seq)
    if rc == seq:
        return False

    return _has_min_hits(rc, rrna_kmers, k, min_hits)
