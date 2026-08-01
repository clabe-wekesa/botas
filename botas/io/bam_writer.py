"""
bam_writer.py

Minimal, SAM-compliant BAM writer for BOTAS.

This module is intentionally simple and tightly coupled to the
current CLI and aligner implementation.

Responsibilities:
- Create a valid SAM/BAM header
- Write mapped and unmapped single-end reads
- Write paired-end reads
- Remain stable and predictable for benchmarking and publication

Important:
- Coordinate normalization for circular padded references must happen
  before calling this writer.
- This writer only writes the coordinates it receives.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re

import pysam


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

logger = logging.getLogger(__name__)


def ref_aligned_length(cigar: str) -> int:
    nref = 0
    for n, op in _CIGAR_RE.findall(cigar):
        if op in ("M", "D", "N", "=", "X"):
            nref += int(n)
    return nref


@dataclass
class BamWriter:
    out: pysam.AlignmentFile
    ref_names: list[str]
    ref_lengths: list[int]
    ref_id: dict[str, int]

    def close(self) -> None:
        self.out.close()


def open_bam_writer(
    out_bam: str,
    *,
    ref_names: list[str],
    ref_lengths: list[int],
) -> BamWriter:
    if len(ref_names) != len(ref_lengths):
        raise ValueError("ref_names and ref_lengths must match")

    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [
            {"SN": n, "LN": int(l)}
            for n, l in zip(ref_names, ref_lengths)
        ],
        "PG": [
            {
                "ID": "botas",
                "PN": "botas",
                "VN": "dev",
                "CL": "seed-and-extend",
            }
        ],
    }

    bf = pysam.AlignmentFile(out_bam, "wb", header=header)
    ref_id = {n: i for i, n in enumerate(ref_names)}

    return BamWriter(
        out=bf,
        ref_names=ref_names,
        ref_lengths=ref_lengths,
        ref_id=ref_id,
    )


def write_hit(
    bw: BamWriter,
    *,
    qname: str,
    seq: str,
    qual: str,
    hit_rname: str,
    pos0: int,
    cigar: str,
    mapq: int,
    strand: str,
    nm: int,
    ascore: int,
    junction: bool = False,
) -> None:
    """
    Write a mapped single-end read.

    Coordinates are assumed to already be normalized to the BAM
    reference coordinate system.
    """
    a = pysam.AlignedSegment(bw.out.header)

    a.query_name = qname
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)

    a.flag = 0
    if strand == "-":
        a.flag |= 0x10

    a.reference_id = bw.ref_id[hit_rname]
    a.reference_start = int(pos0)
    a.mapping_quality = max(0, min(60, int(mapq)))
    a.cigarstring = cigar

    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0

    a.set_tag("NM", int(nm), value_type="i")
    a.set_tag("AS", int(ascore), value_type="i")
    if junction:
        a.set_tag("XC", "JUNCTION", value_type="Z")

    bw.out.write(a)


def write_unmapped(
    bw: BamWriter,
    *,
    qname: str,
    seq: str,
    qual: str,
) -> None:
    """
    Write an unmapped single-end read.
    """
    a = pysam.AlignedSegment(bw.out.header)

    a.query_name = qname
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)

    a.flag = 0x4
    a.reference_id = -1
    a.reference_start = -1
    a.mapping_quality = 0
    a.cigarstring = "*"

    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0

    bw.out.write(a)


def write_pair(
    bw: BamWriter,
    *,
    qname: str,
    r1_seq: str,
    r1_qual: str,
    r2_seq: str,
    r2_qual: str,
    hit1,
    hit2,
    proper_pair: bool,
    insert_size: int | None,
) -> None:
    """
    Write paired-end reads.

    Coordinates are assumed to already be normalized before calling
    this writer.
    """

    def _make(
        *,
        is_read1: bool,
        seq: str,
        qual: str,
        hit,
        mate_hit,
    ) -> pysam.AlignedSegment:
        a = pysam.AlignedSegment(bw.out.header)

        a.query_name = qname
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(qual)

        flag = 0x1
        flag |= 0x40 if is_read1 else 0x80

        if hit is None:
            flag |= 0x4
        elif hit.strand == "-":
            flag |= 0x10

        if mate_hit is None:
            flag |= 0x8
        elif mate_hit.strand == "-":
            flag |= 0x20

        if proper_pair and hit is not None and mate_hit is not None:
            flag |= 0x2

        a.flag = flag
        a.flag &= ~(0x100 | 0x800)

        if hit is None:
            a.reference_id = -1
            a.reference_start = -1
            a.mapping_quality = 0
            a.cigarstring = "*"
        else:
            a.reference_id = bw.ref_id[hit.rname]
            a.reference_start = int(hit.pos0)
            a.mapping_quality = max(0, min(60, int(hit.mapq)))
            a.cigarstring = hit.cigar
            a.set_tag("NM", int(-hit.ascore), value_type="i")
            a.set_tag("AS", int(hit.ascore), value_type="i")
            if getattr(hit, "junction", False):
                a.set_tag("XC", "JUNCTION", value_type="Z")

        if mate_hit is None:
            a.next_reference_id = -1
            a.next_reference_start = -1
        else:
            a.next_reference_id = bw.ref_id[mate_hit.rname]
            a.next_reference_start = int(mate_hit.pos0)

        if hit is not None and mate_hit is not None:
            start1 = int(hit.pos0)
            end1 = start1 + ref_aligned_length(hit.cigar)
            start2 = int(mate_hit.pos0)
            end2 = start2 + ref_aligned_length(mate_hit.cigar)

            left = min(start1, start2)
            right = max(end1, end2)
            tlen = right - left

            if start1 == left:
                a.template_length = int(tlen)
            else:
                a.template_length = -int(tlen)
        else:
            a.template_length = 0

        return a

    a1 = _make(
        is_read1=True,
        seq=r1_seq,
        qual=r1_qual,
        hit=hit1,
        mate_hit=hit2,
    )
    a2 = _make(
        is_read1=False,
        seq=r2_seq,
        qual=r2_qual,
        hit=hit2,
        mate_hit=hit1,
    )

    bw.out.write(a1)
    bw.out.write(a2)
