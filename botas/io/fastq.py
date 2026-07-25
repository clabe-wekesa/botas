# botas/io/fastq.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Generator, Iterable, Optional, TextIO, Tuple
import gzip
import io
import os


@dataclass(frozen=True)
class FastqRead:
    name: str
    seq: str
    plus: str
    qual: str


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fastq(path: str) -> Generator[FastqRead, None, None]:
    """
    Streaming FASTQ reader for .fq/.fastq (optionally gzipped).
    """
    with _open_text(path) as fh:
        while True:
            name = fh.readline().rstrip("\n")
            if not name:
                return
            seq = fh.readline().rstrip("\n")
            plus = fh.readline().rstrip("\n")
            qual = fh.readline().rstrip("\n")
            if not qual:
                return
            if not name.startswith("@"):
                raise ValueError(f"Invalid FASTQ record (name line): {name[:50]}")
            yield FastqRead(name=name[1:].split()[0], seq=seq, plus=plus, qual=qual)
