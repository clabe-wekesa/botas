from botas.operons.features import same_strand, intergenic_distance, coverage_ratio


def same_operon(g1, g2, cov, max_igd, min_coverage, min_cov_ratio=0.5):
    if not same_strand(g1, g2):
        return False

    if intergenic_distance(g1, g2) > max_igd:
        return False

    c1 = cov.get(g1["id"], 0.0)
    c2 = cov.get(g2["id"], 0.0)

    if c1 < min_coverage or c2 < min_coverage:
        return False

    if coverage_ratio(c1, c2) < min_cov_ratio:
        return False

    return True


def operon_score(igds, mean_cov, min_cov, cv, max_igd):
    if not igds or mean_cov <= 0:
        return 0.0
    igd_term = 1.0 - (sum(igds) / len(igds)) / max_igd
    cov_term = min_cov / mean_cov
    stab_term = 1.0 / (1.0 + cv)
    s = igd_term * cov_term * stab_term
    return max(0.0, min(1.0, s))

def operon_confidence(score):
    if score >= 0.8:
        return "HIGH"
    if score >= 0.5:
        return "MEDIUM"
    return "LOW"

