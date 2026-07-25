#!/usr/bin/env python3
"""
rrna.py — rRNA detection helpers for bota
"""

import logging
import gzip
import pysam
from Bio import SeqIO

def load_rrna_kmers(rrna_fa: str, k: int = 17) -> set[str]:
    """
    Build a set of k-mers from an rRNA FASTA file (.fa or .fa.gz).
    Logs the source and total k-mers.
    """
    rrna_kmers = set()
    logging.info(f"[rrna] loading rRNA database → {rrna_fa}")
    try:
        opener = gzip.open if rrna_fa.endswith(".gz") else open
        with opener(rrna_fa, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper().replace("U", "T")
                for i in range(0, max(0, len(seq) - k + 1)):
                    rrna_kmers.add(seq[i:i + k])
    except Exception as e:
        logging.error(f"[rrna] failed to load {rrna_fa}: {e}")
    logging.info(f"[rrna] built k-mer set (k={k}) size={len(rrna_kmers):,}")
    return rrna_kmers


def is_rrna_like(seq: str, rrna_kmers: set[str], k: int = 17, min_hits: int = 3) -> bool:
    """
    Check if a sequence shares ≥min_hits k-mers with the rRNA k-mer set.

    Args:
        seq: input read sequence
        rrna_kmers: set of k-mers derived from rRNA FASTA
        k: k-mer size
        min_hits: minimum number of shared k-mers for rRNA classification

    Returns:
        True if sequence is rRNA-like, False otherwise
    """
    seq = seq.upper()
    hits = 0
    for i in range(len(seq) - k + 1):
        if seq[i:i + k] in rrna_kmers:
            hits += 1
            if hits >= min_hits:
                return True
    return False
