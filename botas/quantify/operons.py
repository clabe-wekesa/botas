# botas/quantify/operons.py
from __future__ import annotations
from typing import Dict, List
import csv
import logging
from typing import Dict, List, Tuple
from botas.quantify.genes import GeneFeature


def compute_operon_metadata(operons: Dict[str, List[str]], genes: Dict[str, GeneFeature]):
    """
    Returns:
      operon_id -> dict with contig, start, end, strand, genes_csv
    Coordinates are 0-based half-open [start, end).
    """
    meta = {}

    for op_id, gids in operons.items():
        gfs = [genes[g] for g in gids if g in genes]

        if not gfs:
            meta[op_id] = {
                "contig": "",
                "start": -1,
                "end": -1,
                "strand": ".",
                "genes": ",".join(gids),
            }
            continue

        contig = gfs[0].contig
        strand = gfs[0].strand

        starts = []
        ends = []

        for gf in gfs:
            for s, e in gf.blocks:
                starts.append(s)
                ends.append(e)

        meta[op_id] = {
            "contig": contig,
            "start": min(starts),
            "end": max(ends),
            "strand": strand,
            "genes": ",".join(gids),
        }

    return meta


def load_operons_tsv(path: str, gene_col: str) -> Dict[str, List[str]]:
    operons: Dict[str, List[str]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []

        if "operon_id" not in cols:
            raise ValueError(f"Operon TSV must contain 'operon_id'. Found: {cols}")

        if gene_col not in cols:
            raise ValueError(
                f"Expected column '{gene_col}' in operon TSV. Found columns: {cols}"
            )

        for row in reader:
            op_id = row["operon_id"]
            raw = row[gene_col]

            if not raw:
                continue

            # 🔴 THIS IS THE CRITICAL FIX
            gids = [g.strip() for g in raw.split(",") if g.strip()]

            operons.setdefault(op_id, []).extend(gids)

    return operons



def aggregate_operon_counts(
    gene_counts: Dict[str, float],
    operons: Dict[str, List[str]],
) -> Dict[str, float]:
    """
    operon_count = sum of member gene counts
    """
    out: Dict[str, float] = {}
    for op_id, genes in operons.items():
        out[op_id] = sum(gene_counts.get(g, 0.0) for g in genes)
    return out


def aggregate_operon_lengths(
    gene_lengths: Dict[str, int],
    operons: Dict[str, List[str]],
) -> Dict[str, int]:
    """
    operon_length = sum of member gene lengths
    """
    out: Dict[str, int] = {}
    for op_id, genes in operons.items():
        out[op_id] = sum(gene_lengths.get(g, 0) for g in genes)
    return out
