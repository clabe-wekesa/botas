from pathlib import Path
import csv

def infer_operons(bam: str, gff: str, min_score: float = 0.5, classifier: str = "rules"):
    """
    Infer operons and return path to operon TSV.
    """

    out_tsv = Path(bam).with_suffix(".operons.tsv")

    print(f"[botas] Writing inferred operons → {out_tsv}")

    # TEMPORARY: one-gene-per-operon fallback
    # (This preserves pipeline correctness while logic evolves)
    genes = []

    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 9:
                continue
            if f[2].lower() != "gene":
                continue

            attrs = f[8]
            gid = None
            if "ID=" in attrs:
                gid = attrs.split("ID=")[1].split(";")[0]
            if gid:
                genes.append(gid)

    with open(out_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["operon_id", "gene_id"])

        for i, gid in enumerate(genes, start=1):
            writer.writerow([f"op{i:05d}", gid])

    return str(out_tsv)
