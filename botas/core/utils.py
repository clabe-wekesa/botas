# botas/core/utils.py
from __future__ import annotations
from typing import Iterable, Iterator, TypeVar

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

T = TypeVar("T")


def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def normalize_read_name(name: str) -> str:
    """
    Return the read name up to the first whitespace.

    FIX: replaces name.split()[0] which scans the entire string and
    allocates a full list. str.partition stops at the first separator
    and allocates only two string objects regardless of name length.
    """
    head, _, _ = name.partition(" ")
    return head


def chunked(iterable: Iterable[T], size: int) -> Iterator[list[T]]:
    """
    Yield successive chunks of size `size` from an iterable.

    Parameters
    ----------
    iterable : Iterable[T]
        Input iterable.
    size : int
        Chunk size (>0).

    Yields
    ------
    list[T]
        Lists of at most `size` elements.
    """
    if size <= 0:
        raise ValueError("chunk size must be > 0")

    chunk: list[T] = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
