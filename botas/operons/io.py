import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from urllib.parse import unquote

import pysam

from botas.operons.classifier import operon_score, operon_confidence


DISPLAY_ID_CHOICES = ("gene", "locus_tag", "Name", "ID")


def _parse_attrs(attr_str):
    """Parse GFF3 attributes into a decoded dictionary."""
    attrs = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def _first_nonempty(attrs, keys):
    for key in keys:
        value = attrs.get(key)
        if value:
            return value
    return None


def _pick_display_id(attrs, preferred="gene"):
    """Choose the user-facing identifier, preferring the requested GFF field."""
    if preferred not in DISPLAY_ID_CHOICES:
        raise ValueError(
            f"Unsupported gene identifier {preferred!r}; "
            f"choose one of {', '.join(DISPLAY_ID_CHOICES)}"
        )

    fallback = ["gene", "locus_tag", "Name", "ID"]
    ordered = [preferred] + [key for key in fallback if key != preferred]
    return _first_nonempty(attrs, ordered) or "unknown"


def _pick_internal_id(attrs, chrom, start, end, strand, row_number):
    """
    Return a stable internal key for coverage and pair bookkeeping.

    This is deliberately independent of the user-facing identifier so that
    duplicate gene symbols (or missing gene names) cannot overwrite coverage.
    """
    base = _first_nonempty(attrs, ("ID", "locus_tag"))
    if base:
        return base
    return f"{chrom}:{start}-{end}:{strand}:row{row_number}"


def load_genes(gff_path, feature_types=("gene",), id_attribute="gene"):
    """
    Load GFF3 gene features.

    ``internal_id`` is used for all computation. ``display_id`` is the value
    shown in output and is controlled by ``id_attribute``.
    """
    genes = []
    seen_internal = set()

    with open(gff_path, "r", encoding="utf-8") as fh:
        for row_number, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, _, ftype, start_s, end_s, _, strand, _, attrs_s = fields
            if ftype not in feature_types:
                continue
            if strand not in {"+", "-"}:
                continue

            start = int(start_s)
            end = int(end_s)
            attrs = _parse_attrs(attrs_s)

            internal_id = _pick_internal_id(
                attrs,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                row_number=row_number,
            )

            # IDs in malformed GFF files are sometimes duplicated. Preserve
            # every feature without allowing one gene to overwrite another.
            if internal_id in seen_internal:
                internal_id = f"{internal_id}#{chrom}:{start}-{end}:{strand}"
            seen_internal.add(internal_id)

            genes.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "internal_id": internal_id,
                    "display_id": _pick_display_id(attrs, preferred=id_attribute),
                    "gene": attrs.get("gene"),
                    "locus_tag": attrs.get("locus_tag"),
                    "gff_id": attrs.get("ID"),
                    "name": attrs.get("Name"),
                }
            )

    return genes


def _cov_one_region(bam_path, contig, region_start, region_end, genes):
    """Calculate partial gene depth within one non-overlapping region."""
    partial_depth = {g["internal_id"]: 0 for g in genes}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if contig not in bam.references:
            return partial_depth

        for col in bam.pileup(
            contig,
            region_start,
            region_end,
            truncate=True,
            stepper="all",
            min_base_quality=0,
        ):
            pos = col.reference_pos + 1
            depth = col.nsegments

            for gene in genes:
                if gene["start"] <= pos <= gene["end"]:
                    partial_depth[gene["internal_id"]] += depth

    return partial_depth


def _make_coverage_tasks(by_contig, window_size):
    tasks = []

    for contig, genes in by_contig.items():
        if not genes:
            continue

        first_position = min(g["start"] for g in genes)
        last_position = max(g["end"] for g in genes)
        region_start = first_position - 1

        while region_start < last_position:
            region_end = min(region_start + window_size, last_position)
            overlapping_genes = [
                g
                for g in genes
                if g["end"] > region_start and g["start"] <= region_end
            ]

            if overlapping_genes:
                tasks.append((contig, region_start, region_end, overlapping_genes))

            region_start = region_end

    return tasks


def compute_gene_coverage(bam_path, genes, max_workers=None, window_size=250_000):
    """Calculate mean per-base coverage for every gene."""
    by_contig = {}
    for gene in genes:
        by_contig.setdefault(gene["chrom"], []).append(gene)

    for contig_genes in by_contig.values():
        contig_genes.sort(key=lambda gene: gene["start"])

    total_depth = {gene["internal_id"]: 0 for gene in genes}
    gene_len = {
        gene["internal_id"]: max(1, gene["end"] - gene["start"] + 1)
        for gene in genes
    }

    tasks = _make_coverage_tasks(by_contig=by_contig, window_size=window_size)

    print(
        f"[coverage] {len(genes)} genes, "
        f"{len(tasks)} regions, "
        f"{max_workers or 'default'} workers"
    )

    if max_workers == 1:
        for contig, start, end, region_genes in tasks:
            partial_depth = _cov_one_region(
                bam_path, contig, start, end, region_genes
            )
            for gene_id, depth in partial_depth.items():
                total_depth[gene_id] += depth
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    _cov_one_region,
                    bam_path,
                    contig,
                    start,
                    end,
                    region_genes,
                )
                for contig, start, end, region_genes in tasks
            ]

            for future in as_completed(futures):
                partial_depth = future.result()
                for gene_id, depth in partial_depth.items():
                    total_depth[gene_id] += depth

    return {
        gene_id: total_depth[gene_id] / gene_len[gene_id]
        for gene_id in total_depth
    }


def operon_stats(genes, cov):
    if not genes:
        return 0.0, 0.0, 0.0

    covs = [cov.get(g["internal_id"], 0.0) for g in genes]
    mean = sum(covs) / len(covs)
    mn = min(covs)

    if mean > 0:
        cv = math.sqrt(sum((c - mean) ** 2 for c in covs) / len(covs)) / mean
    else:
        cv = 0.0

    return mean, mn, cv


def write_operons_gff(path, operons, cov, max_igd, attribute_prefix=""):
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
                f"genes={','.join(g['display_id'] for g in op)};"
                f"{attribute_prefix}mean_cov={mean_cov:.3f};"
                f"{attribute_prefix}score={score:.3f};"
                f"{attribute_prefix}confidence={conf}"
            )

            gff.write(
                f"{chrom}\tbotas\toperon\t{start}\t{end}\t{score:.3f}\t"
                f"{strand}\t.\t{attrs}\n"
            )
