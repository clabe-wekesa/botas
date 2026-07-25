# botas/io/reference_set.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Iterator, Tuple, Optional
import gzip


# ------------------------------------------------------------
# FASTA reader (multi-contig)
# ------------------------------------------------------------

def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    seq_chunks: List[str] = []

    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks).upper()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                if name is None:
                    raise ValueError(f"FASTA parse error in {path}")
                seq_chunks.append(line)

    if name is not None:
        yield name, "".join(seq_chunks).upper()


# ------------------------------------------------------------
# Reference data structures
# ------------------------------------------------------------

@dataclass
class ContigRef:
    name: str
    seq: str
    circular: bool
    index: object | None = None
    orig_len: int | None = None
    circular_overhang: int = 0
    circular_overhang_percent: int = 100

    @property
    def length(self) -> int:
        return len(self.seq)


class ReferenceSet:
    def __init__(self, contigs: List[ContigRef]):
        if not contigs:
            raise ValueError("ReferenceSet: no contigs loaded")

        seen = set()
        for c in contigs:
            if c.name in seen:
                raise ValueError(f"Duplicate contig name: {c.name}")
            seen.add(c.name)

        self._contigs = contigs
        self._by_name: Dict[str, ContigRef] = {c.name: c for c in contigs}

    def contigs(self) -> List[ContigRef]:
        return self._contigs

    def get(self, name: str) -> ContigRef:
        return self._by_name[name]

    def names(self) -> List[str]:
        return [c.name for c in self._contigs]

    def __len__(self) -> int:
        return len(self._contigs)


# ------------------------------------------------------------
# Loader
# ------------------------------------------------------------

def load_reference_set(
    path: str,
    *,
    circular: bool = False,
    circular_contigs: list[str] | None = None,
    circular_overhang_percent: int = 100,
) -> ReferenceSet:

    if not 0 <= int(circular_overhang_percent) <= 100:
        raise ValueError("circular_overhang_percent must be between 0 and 100")
    circular_overhang_percent = int(circular_overhang_percent)

    contigs: List[ContigRef] = []

    for name, seq in iter_fasta(path):
        if not seq:
            raise ValueError(f"Empty contig: {name}")

        circ_set = set(circular_contigs or [])

        if circular:
            is_circ = True
        else:
            is_circ = name in circ_set

        orig_len = len(seq)
        overhang = 0

        if is_circ:
            if circular_overhang_percent == 100:
                # Backward-compatible full circular unrolling.
                seq2 = seq + seq
            elif circular_overhang_percent == 0:
                # Useful as a linear-control setting while preserving metadata.
                seq2 = seq
            else:
                # Partial circular padding: one cheap one-time copy per contig.
                overhang = max(1, int(orig_len * circular_overhang_percent / 100))
                overhang = min(overhang, orig_len)
                seq2 = seq[-overhang:] + seq + seq[:overhang]
        else:
            seq2 = seq

        contigs.append(
            ContigRef(
                name=name,
                seq=seq2,
                circular=is_circ,
                orig_len=orig_len,
                circular_overhang=overhang,
                circular_overhang_percent=circular_overhang_percent if is_circ else 0,
            )
        )

    return ReferenceSet(contigs)
