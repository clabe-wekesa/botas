def merge_pairs(genes, pair_labels):
    if not genes:
        return []

    if len(pair_labels) != len(genes) - 1:
        raise ValueError("pair_labels length must be len(genes) - 1")

    operons = []
    current = [genes[0]]

    for i, keep in enumerate(pair_labels):
        if keep:
            current.append(genes[i + 1])
        else:
            operons.append(current)
            current = [genes[i + 1]]

    operons.append(current)
    return operons
