# botas/core/botas_index.py
from __future__ import annotations

import pickle
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Tuple
from botas.core.utils import revcomp
from botas.io.reference_set import load_reference_set


BOTAS_INDEX_MAGIC = "BOTAS_INDEX_V1"
BOTAS_INDEX_SUFFIX = ".botas.idx"
_BASE_CODE = {"A": 0, "a": 0, "C": 1, "c": 1, "G": 2, "g": 2, "T": 3, "t": 3,}


def default_index_path(ref_fa: str) -> str:
    return str(Path(ref_fa).with_suffix(Path(ref_fa).suffix + BOTAS_INDEX_SUFFIX))


def encode_kmer(seq: str) -> int | None:
    """
    Encode DNA k-mer as 2-bit integer.
    A=0, C=1, G=2, T=3.
    Return None if sequence contains non-ACGT.
    """
    code = 0
    for ch in seq.upper():
        code <<= 2
        if ch == "A":
            code |= 0
        elif ch == "C":
            code |= 1
        elif ch == "G":
            code |= 2
        elif ch == "T":
            code |= 3
        else:
            return None
    return code


def canonical_kmer_code(kmer: str) -> int | None:
    fwd = encode_kmer(kmer)
    if fwd is None:
        return None

    rev = encode_kmer(revcomp(kmer))
    if rev is None:
        return None

    return fwd if fwd < rev else rev


@dataclass
class BotasIndexedContig:
    name: str
    seq: str
    circular: bool
    orig_len: int
    circular_overhang: int
    circular_overhang_percent: int
    k: int
    w: int
    index: Dict[int, array]

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
class BotasIndex:
    ref_fa: str
    k: int
    w: int
    circular_all: bool
    contigs: list[BotasIndexedContig]


class BotasIntIndexAdapter:
    """
    Adapter that behaves like current KmerIndex enough for align_core.py.

    It provides:
      - windows_for_read()
      - iter_read_minimizers()
      - num_seed_hits()

    Internally it uses integer-coded minimizers and compact arrays.
    """

    def __init__(
        self,
        *,
        seq: str,
        index: Dict[int, array],
        k: int,
        w: int,
        circular: bool,
        label: str | None = None,
    ):
        self.ref_seq = seq
        self.ref_len = len(seq)
        self._index = index
        self.k = k
        self.w = w
        self.circular = circular
        self.label = label

    def iter_read_minimizers(
        self,
        read_seq: str,
        step: int = 1,
    ):
        """
        Fast read minimizer generator using integer-encoded k-mers.

        This avoids repeated string slicing + revcomp() calls for every
        k-mer. It computes forward and reverse-complement k-mer codes using
        rolling 2-bit encoding.
        """
        from collections import deque

        k, w = self.k, self.w
        read_len = len(read_seq)

        if read_len < k + w - 1:
            return

        mask = (1 << (2 * k)) - 1

        def base_code(ch: str):
            if ch == "A" or ch == "a":
                return 0
            if ch == "C" or ch == "c":
                return 1
            if ch == "G" or ch == "g":
                return 2
            if ch == "T" or ch == "t":
                return 3
            return None

        total_kmers = read_len - k + 1
        canonical = []

        fwd = 0
        rev = 0
        valid = 0

        for i, ch in enumerate(read_seq):
            b = _BASE_CODE.get(ch)

            if b is None:
                fwd = 0
                rev = 0
                valid = 0
                continue

            fwd = ((fwd << 2) | b) & mask
            rev = (rev >> 2) | ((3 - b) << (2 * (k - 1)))
            valid += 1

            if valid >= k:
                pos = i - k + 1
                code = fwd if fwd < rev else rev
                canonical.append((code, pos))

        if len(canonical) < w:
            return

        num_windows = len(canonical) - w + 1
        dq = deque()
        last_emitted = None
        added_until = -1

        for i in range(0, num_windows, step):
            right = i + w - 1

            while added_until < right:
                added_until += 1
                code, pos = canonical[added_until]

                while dq and canonical[dq[-1]][0] >= code:
                    dq.pop()

                dq.append(added_until)

            while dq and dq[0] < i:
                dq.popleft()

            if not dq:
                continue

            code, pos = canonical[dq[0]]

            item = (pos, code)
            if item == last_emitted:
                continue

            last_emitted = item
            yield item

    def num_seed_hits(self, read_seq: str) -> int:
        total = 0
        for _, code in self.iter_read_minimizers(read_seq):
            positions = self._index.get(code)
            if positions:
                total += len(positions)
        return total

    def windows_for_read(
        self,
        read_seq: str,
        *,
        pad: int = 250,
        k: int | None = None,
        step: int = 5,
        min_hits: int = 2,
        bin_size: int = 20,
        top_bins: int = 6,
        max_positions_per_kmer: int = 50,
    ):
        from botas.core.ref_index import Window, merge_windows

        if k is None:
            k = self.k
        if k != self.k:
            raise ValueError("Read k must match index k")

        read_len = len(read_seq)
        if read_len < k:
            return []

        votes: dict[int, float] = {}
        L = self.ref_len

        for qpos, code in self.iter_read_minimizers(read_seq, step=step):
            positions = self._index.get(code)
            if not positions:
                continue
            if len(positions) > max_positions_per_kmer:
                continue

            for rpos in positions:
                pred_start = rpos - qpos

                if self.circular:
                    pred_start %= L
                    b = pred_start // bin_size
                    votes[b] = votes.get(b, 0) + 1.0

                    if pred_start < pad:
                        votes[(pred_start + L) // bin_size] = votes.get((pred_start + L) // bin_size, 0) + 0.5
                    elif pred_start > L - pad:
                        votes[(pred_start - L) // bin_size] = votes.get((pred_start - L) // bin_size, 0) + 0.5
                else:
                    b = pred_start // bin_size
                    votes[b] = votes.get(b, 0) + 1.0

        if not votes:
            return []

        ranked = sorted(votes.items(), key=lambda x: x[1], reverse=True)

        windows = []
        for b, cnt in ranked[:top_bins]:
            if cnt < min_hits:
                continue

            pred_start = b * bin_size

            if self.circular:
                start0 = pred_start % L
                end0 = start0 + read_len + 2 * pad
                windows.append(Window(start0=start0 - pad, end0=end0, hits=int(cnt)))
            else:
                windows.append(
                    Window(
                        start0=max(0, pred_start - pad),
                        end0=pred_start + read_len + pad,
                        hits=int(cnt),
                    )
                )

        windows = merge_windows(windows, gap=pad)
        return [w for w in windows if w.start0 < w.end0]


def build_botas_index(
    ref_fa: str,
    *,
    k: int = 15,
    w: int = 10,
    circular: bool = False,
    circular_contigs: list[str] | None = None,
    circular_overhang_percent: int = 100,
) -> BotasIndex:
    refset = load_reference_set(
        ref_fa,
        circular=circular,
        circular_contigs=circular_contigs,
        circular_overhang_percent=circular_overhang_percent,
    )

    indexed_contigs = []

    for contig in refset.contigs():
        seq = contig.seq
        L = len(seq)

        index: Dict[int, array] = {}

        if L >= k + w - 1:
            from collections import deque

            limit = L - k + 1
            dq = deque()
            last_added = None

            for pos in range(limit):
                code = canonical_kmer_code(seq[pos:pos + k])

                # Remove invalid or larger kmers from right
                if code is not None:
                    while dq and dq[-1][0] >= code:
                        dq.pop()
                    dq.append((code, pos))

                # Remove kmers outside current minimizer window
                window_start = pos - w + 1
                while dq and dq[0][1] < window_start:
                    dq.popleft()

                # Emit minimizer once the first full window is available
                if window_start >= 0 and dq:
                    best_code, best_pos = dq[0]

                    # Avoid repeated identical minimizer positions
                    if last_added != (best_code, best_pos):
                        if best_code not in index:
                            index[best_code] = array("I")
                        index[best_code].append(best_pos)
                        last_added = (best_code, best_pos)

        indexed_contigs.append(
            BotasIndexedContig(
                name=contig.name,
                seq=contig.seq,
                circular=contig.circular,
                orig_len=contig.orig_len or len(contig.seq),
                circular_overhang=getattr(contig, "circular_overhang", 0),
                circular_overhang_percent=getattr(contig, "circular_overhang_percent", 100 if contig.circular else 0),
                k=k,
                w=w,
                index=index,
            )
        )

    return BotasIndex(
        ref_fa=ref_fa,
        k=k,
        w=w,
        circular_all=circular,
        contigs=indexed_contigs,
    )


def save_botas_index(index: BotasIndex, out_path: str) -> None:
    payload = {
        "magic": BOTAS_INDEX_MAGIC,
        "index": index,
    }

    with open(out_path, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_botas_index(path: str) -> BotasIndex:
    with open(path, "rb") as fh:
        payload = pickle.load(fh)

    if not isinstance(payload, dict) or payload.get("magic") != BOTAS_INDEX_MAGIC:
        raise ValueError(f"Not a valid BOTAS index: {path}")

    return payload["index"]


def attach_adapters(index: BotasIndex) -> list[BotasIndexedContig]:
    """
    Attach adapter indexes to each indexed contig so existing align code
    can use c.index.windows_for_read().
    """
    for c in index.contigs:
        c.index_adapter = BotasIntIndexAdapter(
            seq=c.seq,
            index=c.index,
            k=c.k,
            w=c.w,
            circular=c.circular,
            label=c.name,
        )
        c.index_adapter.orig_len = getattr(c, "orig_len", len(c.seq))
        c.index_adapter.circular_overhang = getattr(c, "circular_overhang", 0)
    return index.contigs