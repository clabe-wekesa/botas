#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
import re


def parse_attrs(text):
    attrs = {}
    for field in text.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k] = unquote(v)
    return attrs


def fasta_lengths(fasta):
    lengths = {}
    name = None
    seq_len = 0

    with open(fasta) as fh:
        for line in fh:
            line = line.strip()

            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len

                name = line[1:].split()[0]
                seq_len = 0

            else:
                seq_len += len(line)

        if name is not None:
            lengths[name] = seq_len

    return lengths


def clean(value, default="-"):
    if not value:
        return default

    value = value.replace("\t", " ").strip()

    return value if value else default


def get_attr(attrs, *names, default="-"):
    for name in names:
        if attrs.get(name):
            return clean(attrs[name])

    return default


def cds_length_aa(start, end):
    nt = abs(end - start) + 1
    return nt // 3


def make_files(sample, fasta, gff, outdir):

    lengths = fasta_lengths(fasta)

    cds = []
    rnas = []

    with open(gff) as fh:

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_text = cols

            start = int(start)
            end = int(end)

            attrs = parse_attrs(attr_text)

            feature_lower = feature.lower()

            # =================================================
            # Protein coding genes -> PTT
            # =================================================

            if feature_lower == "cds":

                gene = get_attr(
                    attrs,
                    "gene",
                    "gene_name",
                    default="-"
                )

                locus = get_attr(
                    attrs,
                    "locus_tag",
                    "ID",
                    "Parent",
                    default="-"
                )

                protein_id = get_attr(
                    attrs,
                    "protein_id",
                    "Name",
                    default="-"
                )

                product = get_attr(
                    attrs,
                    "product",
                    "description",
                    default="-"
                )

                cds.append(
                    (
                        seqid,
                        start,
                        end,
                        strand,
                        cds_length_aa(start, end),
                        protein_id,
                        gene,
                        locus,
                        "-",
                        "-",
                        product,
                    )
                )

            # =================================================
            # Structural RNAs -> RNT
            # =================================================

            elif feature_lower in {
                "rrna",
                "trna",
                "ncrna",
                "tmrna",
                "srna",
                "rna"
            }:

                gene = get_attr(
                    attrs,
                    "gene",
                    "gene_name",
                    default="-"
                )

                locus = get_attr(
                    attrs,
                    "locus_tag",
                    "ID",
                    "Parent",
                    default="-"
                )

                product = get_attr(
                    attrs,
                    "product",
                    "Name",
                    "description",
                    default=feature
                )

                rnas.append(
                    (
                        seqid,
                        start,
                        end,
                        strand,
                        gene,
                        locus,
                        product,
                    )
                )


    # =========================================================
    # Rockhopper expects the legacy NCBI table style.
    #
    # For multi-replicon bacterial genomes, keep everything
    # together. Coordinates correspond to the relevant FASTA
    # contig IDs.
    # =========================================================

    ptt_file = outdir / f"{sample}.ptt"
    rnt_file = outdir / f"{sample}.rnt"


    # =========================================================
    # PTT
    # =========================================================

    with open(ptt_file, "w") as out:

        total_len = sum(lengths.values())

        out.write(
            f"{sample} - 1..{total_len}\n"
        )

        out.write(
            f"{len(cds)} proteins\n"
        )

        out.write(
            "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"
        )

        for (
            seqid,
            start,
            end,
            strand,
            aa_len,
            protein_id,
            gene,
            locus,
            code,
            cog,
            product,
        ) in cds:

            out.write(
                f"{start}..{end}\t"
                f"{strand}\t"
                f"{aa_len}\t"
                f"{protein_id}\t"
                f"{gene}\t"
                f"{locus}\t"
                f"{code}\t"
                f"{cog}\t"
                f"{product}\n"
            )


    # =========================================================
    # RNT
    # =========================================================

    with open(rnt_file, "w") as out:

        total_len = sum(lengths.values())

        out.write(
            f"{sample} - 1..{total_len}\n"
        )

        out.write(
            f"{len(rnas)} RNAs\n"
        )

        out.write(
            "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"
        )

        for (
            seqid,
            start,
            end,
            strand,
            gene,
            locus,
            product,
        ) in rnas:

            length = abs(end - start) + 1

            out.write(
                f"{start}..{end}\t"
                f"{strand}\t"
                f"{length}\t"
                f"-\t"
                f"{gene}\t"
                f"{locus}\t"
                f"-\t"
                f"-\t"
                f"{product}\n"
            )


    print(
        f"{sample}: "
        f"{len(cds)} proteins, "
        f"{len(rnas)} RNAs"
    )

    print(f"  {ptt_file}")
    print(f"  {rnt_file}")


def main():
    import argparse,csv
    ap=argparse.ArgumentParser(description="Create Rockhopper annotation files from FASTA/GFF pairs in a manifest.")
    ap.add_argument("--manifest",required=True,help="TSV columns: sample,reference,gff")
    ap.add_argument("--outdir",required=True)
    args=ap.parse_args()
    outroot=Path(args.outdir); outroot.mkdir(parents=True,exist_ok=True)
    with open(args.manifest,newline="") as fh:
        r=csv.DictReader(fh,delimiter="\t")
        req={"sample","reference","gff"}
        if not req.issubset(r.fieldnames or []): ap.error("manifest must contain sample, reference, gff")
        for row in r:
            sample=row["sample"].strip()
            if not sample: continue
            fasta=Path(row["reference"]).expanduser(); gff=Path(row["gff"]).expanduser()
            if not fasta.is_file() or not gff.is_file(): raise FileNotFoundError(f"{sample}: missing reference or GFF")
            outdir=outroot/sample; outdir.mkdir(parents=True,exist_ok=True)
            dest_fasta=outdir/f"{sample}.fna"
            if not dest_fasta.exists(): dest_fasta.write_bytes(fasta.read_bytes())
            make_files(sample,fasta,gff,outdir)

if __name__ == "__main__":
    main()
