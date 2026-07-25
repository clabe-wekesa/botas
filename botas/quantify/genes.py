# botas/quantify/genes.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneFeature:
    gene_id: str
    contig: str
    strand: str  # "+", "-"
    blocks: Tuple[Tuple[int, int], ...]  # 0-based, half-open [start,end)
    length: int  # merged exonic length


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def _pick_gene_id(attrs: str, key: str) -> Optional[str]:
    # Robust-enough GFF3 attribute parsing for common cases:
    # ID=gene123;Name=...  or  gene_id=xxx;  or  Parent=...
    # We prefer ID=... then gene_id=...
    parts = attrs.split(";")
    kv = {}
    for p in parts:
        p = p.strip()
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        kv[k.strip()] = v.strip()

    if key == "auto":
        return kv.get("Parent") or kv.get("ID") or kv.get("gene_id")

    return kv.get(key)


def load_gene_features_from_gff(
    gff_path: str,
    feature_type: str = "CDS",
    id_from: str = "auto",
) -> Dict[str, GeneFeature]:
    """
    Phase-1 assumption (bacterial): quantify on CDS blocks grouped by gene.

    Strategy:
      - group intervals by gene_id
      - merge overlapping/adjacent intervals => exon/CDS blocks
      - compute merged length

    Notes:
      - GFF coordinates are 1-based inclusive. We convert to 0-based half-open.
      - For bacterial genomes, CDS entries often correspond to genes (or multiple CDS per gene for split annotations).
    """
    by_gene: Dict[str, List[Tuple[str, str, int, int]]] = {}  # gene_id -> list(contig,strand,s,e)

    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            contig, _source, ftype, start, end, _score, strand, _phase, attrs = cols
            if ftype != feature_type:
                continue

            gid = _pick_gene_id(attrs, id_from)
            if gid is None:
                # fallback: try Parent=... if ID absent (common when CDS has Parent=gene)
                # We won't over-engineer; keep it simple.
                parts = attrs.split(";")
                parent = None
                for p in parts:
                    p = p.strip()
                    if p.startswith("Parent="):
                        parent = p.split("=", 1)[1].strip()
                        break
                gid = parent

            if gid is None:
                continue

            s1 = int(start)
            e1 = int(end)
            # convert 1-based inclusive to 0-based half-open
            s0 = s1 - 1
            e0 = e1  # inclusive -> half-open end = end
            by_gene.setdefault(gid, []).append((contig, strand, s0, e0))

    features: Dict[str, GeneFeature] = {}
    for gid, rows in by_gene.items():
        # genes should not span contigs/strands; keep first, warn if inconsistent
        contig0, strand0, _, _ = rows[0]
        intervals = []
        for contig, strand, s, e in rows:
            if contig != contig0 or strand != strand0:
                logger.warning("Gene %s has inconsistent contig/strand in GFF; keeping first (%s,%s)", gid, contig0, strand0)
                continue
            intervals.append((s, e))
        merged = _merge_intervals(intervals)
        length = sum(e - s for s, e in merged)
        features[gid] = GeneFeature(
            gene_id=gid,
            contig=contig0,
            strand=strand0,
            blocks=tuple(merged),
            length=length,
        )

    logger.info("Loaded %d gene features from %s (type=%s)", len(features), gff_path, feature_type)
    return features
