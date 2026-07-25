from __future__ import annotations
from collections import defaultdict, deque
from botas.core.utils import revcomp
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterator
import multiprocessing as mp
import logging

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

logger = logging.getLogger(__name__)

# NOTE: _RC was a duplicate of utils._COMP and was never used here — removed.


@dataclass(frozen=True)
class Window:
    start0: int
    end0: int
    hits: int = 1


def merge_windows(windows: List[Window], gap: int = 50) -> List[Window]:
    if not windows:
        return []
    windows = sorted(windows, key=lambda w: (w.start0, w.end0))
    merged = [windows[0]]
    for w in windows[1:]:
        last = merged[-1]
        if w.start0 <= last.end0 + gap:
            merged[-1] = Window(
                start0=min(last.start0, w.start0),
                end0=max(last.end0, w.end0),
                hits=last.hits + w.hits,
            )
        else:
            merged.append(w)
    return merged


# ------------------------------------------------------------------
# KmerIndex
# ------------------------------------------------------------------

class KmerIndex:
    """
    Canonical minimizer-based reference index.
    Supports circular genomes.
    """

    def __init__(
        self,
        ref_seq: str,
        *,
        k: int = 15,
        w: int = 10,
        circular: bool = False,
        label: str | None = None,
    ):
        self.k = k
        self.w = w
        self.circular = circular
        self.ref_seq = ref_seq
        self.ref_len = len(ref_seq)
        self._index: Dict[str, List[int]] = {}
        self.label = label
        self._build()

    # ------------------------------------------------------------------
    # Build minimizer index
    # FIX: replaced O(w) inner scan with a deque-based sliding-window
    # minimum, reducing index build from O(L*w) to O(L).
    # FIX: replaced dict.setdefault([], ...) with defaultdict(list) to
    # avoid constructing a throw-away list object on every existing key.
    # ------------------------------------------------------------------

    def _build(self) -> None:
        k, w = self.k, self.w
        seq = self.ref_seq
        L = len(seq)

        if L < k:
            logger.warning("Reference too short for k=%d", k)
            return

        if self.circular:
            seq = seq + seq[:k + w - 1]
            limit = L
        else:
            limit = L - k + 1

        index: Dict[str, List[int]] = defaultdict(list)

        # --- precompute canonical k-mers for every position in seq ---
        total_kmers = len(seq) - k + 1
        canonical: list[tuple[str, int]] = []  # (canonical_kmer, original_pos)
        for i in range(total_kmers):
            kf = seq[i:i + k]
            if "N" in kf:
                canonical.append(("", i))
                continue
            kr = revcomp(kf)
            canonical.append((kf if kf < kr else kr, i))

        # --- sliding-window minimum over the canonical array ---
        # For each window position i in [0, limit), the minimizer is
        # min(canonical[i], ..., canonical[i+w-1]).
        # We use a monotonic deque (indices into `canonical`), keeping
        # only positions where the canonical kmer is non-empty and ≤ all
        # later entries still in the window.

        dq: deque[int] = deque()  # indices into `canonical`, front = current min

        it = range(limit)

        is_main = mp.current_process().name == "MainProcess"
        if (
            tqdm is not None
            and is_main
            and logging.getLogger().isEnabledFor(logging.INFO)
        ):
            desc = "Building minimizer index"
            if self.label:
                desc += f" [{self.label}]"
            it = tqdm(
                it,
                total=limit,
                desc=desc,
                unit="win",
                leave=False,
                mininterval=0.5,
                smoothing=0.0,
                dynamic_ncols=True,
            )

        for i in it:
            # Add the new kmer entering the right side of the window
            new_pos = i + w - 1
            if new_pos < total_kmers:
                new_kmer, _ = canonical[new_pos]
                if new_kmer:  # skip N-containing positions
                    # Pop from back while back is worse (≥) than new entry
                    while dq and (canonical[dq[-1]][0] == "" or canonical[dq[-1]][0] >= new_kmer):
                        dq.pop()
                    dq.append(new_pos)

            # Remove indices that have slid out of the window
            while dq and dq[0] < i:
                dq.popleft()

            if not dq:
                continue

            min_kmer, min_pos = canonical[dq[0]]
            if not min_kmer:
                continue

            pos = min_pos % L if self.circular else min_pos
            index[min_kmer].append(pos)

        self._index = dict(index)
        logger.debug(
            "Minimizer index built: k=%d w=%d circular=%s unique_kmers=%d total_entries=%d",
            k, w, self.circular,
            len(self._index),
            sum(len(v) for v in self._index.values()),
        )

    # ------------------------------------------------------------------
    # Window projection from read
    # FIX: same deque-based O(1) sliding minimizer replaces O(w) inner loop.
    # ------------------------------------------------------------------

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
    ) -> List[Window]:

        if k is None:
            k = self.k
        if k != self.k:
            raise ValueError("Read k must match index k")

        read_len = len(read_seq)
        if read_len < k:
            return []

        votes: Dict[int, float] = {}
        L = self.ref_len

        for qpos, kmer in self.iter_read_minimizers(read_seq, step=step):
            positions = self._index.get(kmer)
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

        windows: List[Window] = []
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
        windows = [w for w in windows if w.start0 < w.end0]

        return windows

    # ------------------------------------------------------------------
    # Read minimizers
    # FIX: deque-based O(1) sliding-window minimum replaces O(w) inner loop.
    # ------------------------------------------------------------------

    def iter_read_minimizers(
        self,
        read_seq: str,
        step: int = 1,
    ) -> Iterator[Tuple[int, str]]:

        k, w = self.k, self.w
        read_len = len(read_seq)

        if read_len < k + w - 1:
            return

        # Precompute canonical k-mers for every position in read_seq
        total_kmers = read_len - k + 1
        canonical: list[tuple[str, int]] = []
        for i in range(total_kmers):
            kf = read_seq[i:i + k]
            if "N" in kf:
                canonical.append(("", i))
                continue
            kr = revcomp(kf)
            canonical.append((kf if kf < kr else kr, i))

        num_windows = read_len - (w + k - 1) + 1
        if num_windows <= 0:
            return

        dq: deque[int] = deque()
        last_emitted: tuple[int, str] | None = None

        for i in range(0, num_windows, step):
            # Add new kmer(s) entering the window on the right
            # When step > 1 we may need to add multiple entries
            new_right = i + w - 1
            # Ensure all positions up to new_right are considered
            # (previous iteration covered up to (i - step) + w - 1)
            first_new = (i - step + w) if i > 0 else 0
            for new_pos in range(max(first_new, 0), min(new_right + 1, total_kmers)):
                new_kmer, _ = canonical[new_pos]
                if new_kmer:
                    while dq and (canonical[dq[-1]][0] == "" or canonical[dq[-1]][0] >= new_kmer):
                        dq.pop()
                    dq.append(new_pos)

            # Expire entries that left the window
            while dq and dq[0] < i:
                dq.popleft()

            if not dq:
                continue

            min_kmer, min_pos = canonical[dq[0]]
            if not min_kmer:
                continue

            # Deduplicate: skip if same (pos, kmer) as last emission
            if last_emitted == (min_pos, min_kmer):
                continue

            last_emitted = (min_pos, min_kmer)
            yield min_pos, min_kmer

    # ------------------------------------------------------------------
    # Seed-hit count helper (used by pe_pool for contig pre-filtering)
    # ------------------------------------------------------------------

    def num_seed_hits(self, read_seq: str) -> int:
        """Count total index hits for a read (fast pre-filter)."""
        total = 0
        for _, kmer in self.iter_read_minimizers(read_seq):
            positions = self._index.get(kmer)
            if positions:
                total += len(positions)
        return total
