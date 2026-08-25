from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

# bam_writer only needs these attributes while the pure TLEN helper is tested.
# The real runtime supplies pysam.
if "pysam" not in sys.modules:
    fake_pysam = types.ModuleType("pysam")
    fake_pysam.AlignmentFile = object
    fake_pysam.AlignedSegment = object
    fake_pysam.qualitystring_to_array = lambda value: value
    sys.modules["pysam"] = fake_pysam

from botas.core.align_core import Hit
from botas.core.align_pe import align_pair_simple
from botas.core.pairing import is_proper_pair_unified
from botas.io.bam_writer import pair_template_lengths


def _hit(pos0: int, strand: str, *, score: int = 0) -> Hit:
    return Hit(
        rname="chr",
        pos0=pos0,
        strand=strand,
        cigar="100M",
        ascore=score,
        mapq=60,
    )


class PairOrientationTests(unittest.TestCase):
    def test_r2_is_aligned_in_its_original_orientation(self) -> None:
        r1 = "A" * 100
        r2 = "C" * 100

        with patch(
            "botas.core.align_pe.align_read",
            side_effect=[_hit(400, "-"), _hit(100, "+")],
        ) as align_read:
            result = align_pair_simple(
                r1_seq=r1,
                r2_seq=r2,
                rname="chr",
                ref_seq="A" * 1000,
                index=object(),
                circular=False,
                max_insert=1200,
                do_rescue=False,
            )

        self.assertEqual(align_read.call_args_list[1].kwargs["read_seq"], r2)
        self.assertEqual(result.hit2.strand, "+")
        self.assertTrue(result.proper_pair)
        self.assertEqual(result.insert_size, 400)

    def test_rescued_mate_strand_depends_on_anchor_strand(self) -> None:
        cases = (
            ([_hit(400, "-"), None], _hit(100, "+"), "+"),
            ([None, _hit(100, "+")], _hit(400, "-"), "-"),
        )

        for initial_hits, rescued_hit, expected_strand in cases:
            with self.subTest(anchor=initial_hits):
                with patch(
                    "botas.core.align_pe.align_read",
                    side_effect=initial_hits,
                ), patch(
                    "botas.core.align_pe.rescue_mate",
                    return_value=rescued_hit,
                ) as rescue:
                    result = align_pair_simple(
                        r1_seq="A" * 100,
                        r2_seq="C" * 100,
                        rname="chr",
                        ref_seq="A" * 1000,
                        index=object(),
                        circular=False,
                        max_insert=1200,
                        do_rescue=True,
                    )

                self.assertEqual(rescue.call_args.kwargs["mate_strand"], expected_strand)
                self.assertTrue(result.proper_pair)


class CircularPairingTests(unittest.TestCase):
    def test_expected_insert_does_not_reject_valid_circular_pair(self) -> None:
        proper, insert_size, orientation = is_proper_pair_unified(
            pos1=900,
            strand1="+",
            pos2=300,
            strand2="-",
            read_len=100,
            circular=True,
            ref_len=1000,
            max_insert=1200,
            expected_insert=300,
        )

        self.assertTrue(proper)
        self.assertEqual(insert_size, 500)
        self.assertEqual(orientation, "FR")

    def test_max_insert_rejects_oversized_circular_pair(self) -> None:
        proper, insert_size, orientation = is_proper_pair_unified(
            pos1=900,
            strand1="+",
            pos2=300,
            strand2="-",
            read_len=100,
            circular=True,
            ref_len=1000,
            max_insert=450,
            expected_insert=300,
        )

        self.assertFalse(proper)
        self.assertEqual(insert_size, 500)
        self.assertEqual(orientation, "FR")


class TemplateLengthTests(unittest.TestCase):
    def test_origin_crossing_pair_uses_circular_insert_size(self) -> None:
        hit1 = SimpleNamespace(rname="chr", pos0=950, cigar="100M")
        hit2 = SimpleNamespace(rname="chr", pos0=50, cigar="100M")

        self.assertEqual(
            pair_template_lengths(hit1=hit1, hit2=hit2, insert_size=200),
            (-200, 200),
        )

    def test_linear_pair_falls_back_to_outer_span(self) -> None:
        hit1 = SimpleNamespace(rname="chr", pos0=100, cigar="90M10S")
        hit2 = SimpleNamespace(rname="chr", pos0=400, cigar="100M")

        self.assertEqual(
            pair_template_lengths(hit1=hit1, hit2=hit2, insert_size=None),
            (400, -400),
        )


if __name__ == "__main__":
    unittest.main()
