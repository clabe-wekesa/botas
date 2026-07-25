from botas.operons.classifier import same_operon
from botas.operons.merge import merge_pairs


def pair_key(g1, g2):
    return (g1["id"], g2["id"])


def build_pair_support(directons, cov_by_bam, args):
    pair_counts = {}
    n_bams = len(cov_by_bam)

    for cov in cov_by_bam:
        for ds in directons:
            for g1, g2 in zip(ds[:-1], ds[1:]):
                if same_operon(
                    g1, g2, cov,
                    max_igd=args.max_igd,
                    min_coverage=args.min_coverage,
                    min_cov_ratio=args.min_cov_ratio,
                ):
                    k = pair_key(g1, g2)
                    pair_counts[k] = pair_counts.get(k, 0) + 1

    return {k: v / n_bams for k, v in pair_counts.items()}


def merge_by_support(directon, pair_support, min_support):
    labels = []
    for g1, g2 in zip(directon[:-1], directon[1:]):
        s = pair_support.get((g1["id"], g2["id"]), 0.0)
        labels.append(s >= min_support)
    return merge_pairs(directon, labels)