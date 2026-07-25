import pysam
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from botas.operons.classifier import operon_score, operon_confidence


def _parse_attrs(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        d[k.strip()] = v.strip()
    return d


def _pick_id(attrs):
    return (
        attrs.get("ID")
        or attrs.get("locus_tag")
        or attrs.get("gene")
        or attrs.get("Name")
        or "unknown"
    )


def load_genes(gff_path, feature_types=("gene",)):
    genes = []
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs_s = fields
            if ftype not in feature_types:
                continue
            attrs = _parse_attrs(attrs_s)
            gid = _pick_id(attrs)
            genes.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "id": gid,
                }
            )
    return genes


def _cov_one_contig(bam_path, contig, gs):
    total_depth = {g["id"]: 0 for g in gs}
    gene_len = {g["id"]: max(1, g["end"] - g["start"] + 1) for g in gs}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if contig not in bam.references:
            return total_depth, gene_len

        print(f"[coverage] {contig}: {len(gs)} genes")

        i = 0
        n = len(gs)
        for col in bam.pileup(
            contig,
            truncate=True,
            stepper="all",
            min_base_quality=0,
        ):
            pos = col.reference_pos + 1
            depth = col.nsegments

            while i < n and gs[i]["end"] < pos:
                i += 1

            j = i
            while j < n and gs[j]["start"] <= pos <= gs[j]["end"]:
                total_depth[gs[j]["id"]] += depth
                j += 1

    return total_depth, gene_len


def compute_gene_coverage(bam_path, genes, max_workers=None):
    # group & sort
    by_contig = {}
    for g in genes:
        by_contig.setdefault(g["chrom"], []).append(g)
    for gs in by_contig.values():
        gs.sort(key=lambda x: x["start"])

    total_depth = {}
    gene_len = {}

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(_cov_one_contig, bam_path, contig, gs): contig
            for contig, gs in by_contig.items()
        }
        for fut in as_completed(futures):
            td, gl = fut.result()
            total_depth.update(td)
            gene_len.update(gl)

    return {gid: total_depth[gid] / gene_len[gid] for gid in total_depth}


def operon_stats(genes, cov):
    if not genes:
        return 0.0, 0.0, 0.0

    covs = [cov.get(g["id"], 0.0) for g in genes]
    mean = sum(covs) / len(covs)
    mn = min(covs)

    if mean > 0:
        cv = math.sqrt(sum((c - mean) ** 2 for c in covs) / len(covs)) / mean
    else:
        cv = 0.0

    return mean, mn, cv


def write_operons_gff(path, operons, cov, max_igd):
    from botas.operons.features import operon_igds
    from botas.operons.classifier import operon_score

    with open(path, "w", encoding="utf-8") as gff:
        gff.write("##gff-version 3\n")
        for i, op in enumerate(operons, 1):
            chrom = op[0]["chrom"]
            strand = op[0]["strand"]
            start = min(g["start"] for g in op)
            end = max(g["end"] for g in op)

            igds = operon_igds(op)
            mean_cov, min_cov, cv = operon_stats(op, cov)
            score = operon_score(igds, mean_cov, min_cov, cv, max_igd)
            conf = operon_confidence(score)

            attrs = (
                f"ID=operon_{i};"
                f"n_genes={len(op)};"
                f"genes={','.join(g['id'] for g in op)};"
                f"mean_cov={mean_cov:.3f};"
                f"score={score:.3f};"
                f"confidence={conf}"
            )

            gff.write(
                f"{chrom}\tbotas\toperon\t{start}\t{end}\t{score:.3f}\t"
                f"{strand}\t.\t{attrs}\n"
            )
