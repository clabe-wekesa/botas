#!/usr/bin/env python3

import csv
import html
import re
from pathlib import Path


# ============================================================
# INPUTS
# ============================================================
import argparse
_parser = argparse.ArgumentParser(description="Build an experimentally supported RegulonDB operon truth set against a supplied GFF3 annotation.")
_parser.add_argument("--operonset", required=True, help="RegulonDB OperonSet TSV")
_parser.add_argument("--gff", required=True, help="GFF3 annotation used for the benchmark")
_parser.add_argument("--operons-out", required=True)
_parser.add_argument("--pairs-out", required=True)
_parser.add_argument("--unmapped-out", required=True)
_args = _parser.parse_args()
OPERONSET = Path(_args.operonset)
GFF = Path(_args.gff)
OUT_OPERONS = Path(_args.operons_out)
OUT_PAIRS = Path(_args.pairs_out)
OUT_UNMAPPED = Path(_args.unmapped_out)


# ============================================================
# HELPERS
# ============================================================

def parse_attrs(text):
    attrs = {}

    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value

    return attrs


def clean_name(name):
    """
    Normalize simple HTML entities/tags that may appear in RegulonDB.
    """

    name = html.unescape(name.strip())

    # remove HTML tags
    name = re.sub(r"<[^>]+>", "", name)

    return name.strip()


# ============================================================
# BUILD CURRENT E. COLI PROTEIN-CODING GENE NAME SET
# ============================================================

gene_id_to_name = {}
protein_coding_ids = set()


with open(GFF) as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        feature = fields[2]
        attrs = parse_attrs(fields[8])

        if feature == "gene":

            gid = attrs.get("ID")

            name = (
                attrs.get("gene")
                or attrs.get("Name")
                or attrs.get("locus_tag")
                or gid
            )

            if gid and name:
                gene_id_to_name[gid] = name

        elif feature == "CDS":

            parent = attrs.get("Parent")

            if parent:
                for p in parent.split(","):
                    protein_coding_ids.add(p)


protein_coding_names = {
    gene_id_to_name[x]
    for x in protein_coding_ids
    if x in gene_id_to_name
}


print(
    f"Current E. coli protein-coding genes: "
    f"{len(protein_coding_names)}"
)


# ============================================================
# FIND REAL OPERONSET HEADER
# ============================================================

with open(
    OPERONSET,
    encoding="utf-8"
) as fh:

    lines = fh.readlines()


header_index = None

for i, line in enumerate(lines):

    if line.startswith("1)operonId\t"):
        header_index = i
        break


if header_index is None:
    raise RuntimeError(
        "Could not locate OperonSet.tsv header."
    )


reader = csv.DictReader(
    lines[header_index:],
    delimiter="\t"
)


# ============================================================
# PROCESS OPERONS
# ============================================================

clean_operons = []
truth_pairs = {}
unmapped_members = []


for row in reader:

    operon_id = row["1)operonId"].strip()
    operon_name = row["2)operonName"].strip()

    n_genes_text = row["6)numberOfGenes"].strip()
    genes_text = row["7)operonGenes"].strip()
    evidence = row["8)operonEvidence"].strip()
    confidence = row["9)confidenceLevel "].strip()


    if not n_genes_text:
        continue

    try:
        n_genes = int(n_genes_text)
    except ValueError:
        continue


    # ========================================================
    # Multi-gene operons only
    # ========================================================

    if n_genes < 2:
        continue


    # ========================================================
    # Experimental evidence only
    # ========================================================

    if "EXP-" not in evidence:
        continue


    raw_genes = [
        clean_name(x)
        for x in genes_text.split("|")
        if x.strip()
    ]


    if len(raw_genes) < 2:
        continue


    # ========================================================
    # Record which members map to current protein-coding genes
    # ========================================================

    valid_mask = [
        gene in protein_coding_names
        for gene in raw_genes
    ]


    for gene, valid in zip(
        raw_genes,
        valid_mask
    ):

        if not valid:

            unmapped_members.append({
                "operon_id": operon_id,
                "operon_name": operon_name,
                "member": gene,
                "evidence": evidence,
            })


    # ========================================================
    # IMPORTANT:
    #
    # Preserve ORIGINAL adjacency.
    #
    # We do NOT first remove invalid members and reconnect
    # remaining genes.
    #
    # Example:
    #
    # geneA | ncRNA | geneB
    #
    # must NOT become:
    #
    # geneA - geneB
    # ========================================================

    valid_pairs_this_operon = []

    for i in range(len(raw_genes) - 1):

        gene1 = raw_genes[i]
        gene2 = raw_genes[i + 1]

        if (
            gene1 not in protein_coding_names
            or gene2 not in protein_coding_names
        ):
            continue

        valid_pairs_this_operon.append(
            (gene1, gene2)
        )


    if not valid_pairs_this_operon:
        continue


    # ========================================================
    # Save cleaned operon information
    # ========================================================

    clean_operons.append({
        "operon_id": operon_id,
        "operon_name": operon_name,
        "original_members": "|".join(raw_genes),
        "n_original_members": len(raw_genes),
        "n_valid_pairs": len(valid_pairs_this_operon),
        "evidence": evidence,
        "confidence": confidence,
    })


    # ========================================================
    # Add truth pairs
    # ========================================================

    for gene1, gene2 in valid_pairs_this_operon:

        # Treat same-operon adjacency as undirected
        key = tuple(sorted((gene1, gene2)))

        if key not in truth_pairs:

            truth_pairs[key] = {
                "gene1": gene1,
                "gene2": gene2,
                "operons": set(),
                "operon_ids": set(),
                "evidence": set(),
                "confidence": set(),
            }


        truth_pairs[key]["operons"].add(
            operon_name
        )

        truth_pairs[key]["operon_ids"].add(
            operon_id
        )

        truth_pairs[key]["evidence"].add(
            evidence
        )

        if confidence:
            truth_pairs[key]["confidence"].add(
                confidence
            )


# ============================================================
# WRITE CLEAN OPERONS
# ============================================================

with open(
    OUT_OPERONS,
    "w",
    newline="",
    encoding="utf-8"
) as out:

    fields = [
        "operon_id",
        "operon_name",
        "original_members",
        "n_original_members",
        "n_valid_pairs",
        "evidence",
        "confidence",
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(clean_operons)


# ============================================================
# WRITE CLEAN TRUTH PAIRS
# ============================================================

with open(
    OUT_PAIRS,
    "w",
    newline="",
    encoding="utf-8"
) as out:

    writer = csv.writer(
        out,
        delimiter="\t"
    )

    writer.writerow([
        "gene1",
        "gene2",
        "operon",
        "operon_id",
        "evidence",
        "confidence",
    ])


    for key in sorted(truth_pairs):

        x = truth_pairs[key]

        writer.writerow([
            x["gene1"],
            x["gene2"],
            ";".join(
                sorted(x["operons"])
            ),
            ";".join(
                sorted(x["operon_ids"])
            ),
            "; ".join(
                sorted(x["evidence"])
            ),
            ";".join(
                sorted(x["confidence"])
            ),
        ])


# ============================================================
# WRITE UNMAPPED MEMBERS
# ============================================================

with open(
    OUT_UNMAPPED,
    "w",
    newline="",
    encoding="utf-8"
) as out:

    fields = [
        "operon_id",
        "operon_name",
        "member",
        "evidence",
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(unmapped_members)


# ============================================================
# SUMMARY
# ============================================================

unique_unmapped = {
    row["member"]
    for row in unmapped_members
}


print()
print(
    "========================================"
)
print(
    "REGULONDB E. COLI TRUTH SUMMARY"
)
print(
    "========================================"
)

print(
    f"Experimental operons with usable pairs: "
    f"{len(clean_operons)}"
)

print(
    f"Protein-coding adjacent truth pairs:     "
    f"{len(truth_pairs)}"
)

print(
    f"Unique unmapped/special members:         "
    f"{len(unique_unmapped)}"
)

print()
print(f"Written: {OUT_OPERONS}")
print(f"Written: {OUT_PAIRS}")
print(f"Written: {OUT_UNMAPPED}")
