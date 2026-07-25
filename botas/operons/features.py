def same_strand(g1, g2):
    return g1["strand"] == g2["strand"]


def intergenic_distance(g1, g2):
    s = g1["strand"]
    if s == "+":
        return g2["start"] - g1["end"] - 1
    if s == "-":
        return g1["start"] - g2["end"] - 1
    raise ValueError(f"Invalid strand: {s!r}")


def coverage_ratio(c1, c2):
    m = max(c1, c2)
    if m <= 0:
        return 0.0
    return min(c1, c2) / m


def operon_igds(genes):
    return [
        intergenic_distance(g1, g2)
        for g1, g2 in zip(genes[:-1], genes[1:])
    ]
