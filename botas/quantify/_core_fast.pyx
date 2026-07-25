# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: nonecheck=False

import pysam


cdef inline bint _strand_matches(
    object aln,
    object gene_strand,
    int strand_mode,
):
    if strand_mode == 0:
        return True

    if strand_mode == 1:
        if aln.is_reverse:
            return gene_strand == "-"
        return gene_strand == "+"

    if strand_mode == 2:
        if aln.is_reverse:
            return gene_strand != "-"
        return gene_strand != "+"

    return False


cdef inline Py_ssize_t _bisect_left(
    tuple starts,
    int value,
):
    cdef Py_ssize_t low = 0
    cdef Py_ssize_t high = len(starts)
    cdef Py_ssize_t middle

    while low < high:
        middle = (low + high) >> 1

        if starts[middle] < value:
            low = middle + 1
        else:
            high = middle

    return low


cdef bint _genes_overlapping_block_fast(
    int block_start,
    int block_end,
    object contig_index,
    object aln,
    int strand_mode,
    set hits,
):
    cdef tuple intervals
    cdef tuple starts
    cdef tuple prefix_max_ends
    cdef tuple interval

    cdef Py_ssize_t position
    cdef int start
    cdef int end
    cdef object gene_id
    cdef object gene_strand

    if block_end <= block_start:
        return False

    intervals = contig_index.intervals
    starts = contig_index.starts
    prefix_max_ends = contig_index.prefix_max_ends

    position = _bisect_left(starts, block_end) - 1

    while position >= 0:
        if prefix_max_ends[position] <= block_start:
            break

        interval = intervals[position]

        start = interval[0]
        end = interval[1]
        gene_id = interval[2]
        gene_strand = interval[3]

        if end > block_start and start < block_end:
            if _strand_matches(
                aln,
                gene_strand,
                strand_mode,
            ):
                hits.add(gene_id)

                if len(hits) > 1:
                    return True

        position -= 1

    return False


cdef tuple _assign_alignment_to_gene_fast(
    object aln,
    object contig_index,
    int strand_mode,
):
    cdef set hits = set()
    cdef object blocks
    cdef object block

    cdef int block_start
    cdef int block_end

    blocks = aln.get_blocks()

    if not blocks:
        return None, False

    for block in blocks:
        block_start = block[0]
        block_end = block[1]

        if _genes_overlapping_block_fast(
            block_start,
            block_end,
            contig_index,
            aln,
            strand_mode,
            hits,
        ):
            return None, True

    if len(hits) == 1:
        return next(iter(hits)), False

    return None, False


cdef inline double _read_weight_fast(
    object aln,
    int multi_mode,
    bint use_nh_tag,
):
    """
    multi_mode:
        0 = ignore multimappers
        1 = unique only
        2 = fractional
    """

    cdef object value
    cdef int nh

    if not use_nh_tag:
        return 1.0

    if not aln.has_tag("NH"):
        return 1.0

    value = aln.get_tag("NH")

    try:
        nh = int(value)
    except (TypeError, ValueError):
        return 1.0

    if multi_mode == 0:
        if nh > 1:
            return 0.0
        return 1.0

    if multi_mode == 1:
        if nh == 1:
            return 1.0
        return 0.0

    if multi_mode == 2:
        if nh <= 1:
            return 1.0
        return 1.0 / nh

    return 1.0


cpdef tuple quantify_bam_fast(
    object bam_path,
    object gene_ids,
    object gene_index,
    int strand_mode,
    int mapq_min,
    int multi_mode,
    bint use_nh_tag,
):
    """
    Quantify one BAM using one compiled execution path.

    Returns
    -------
    counts
        Dictionary mapping gene IDs to weighted counts.

    raw_stats
        Tuple containing:
        total_records,
        primary_mapped_records,
        below_mapq,
        multimapping_discarded,
        no_annotated_contig,
        no_feature,
        ambiguous,
        assigned_records,
        assigned_weight
    """

    cdef dict counts = dict.fromkeys(gene_ids, 0.0)
    cdef dict index_by_tid = {}

    cdef object bam
    cdef object aln
    cdef object contig
    cdef object contig_index
    cdef object gene_id
    cdef object result

    cdef int tid
    cdef double weight
    cdef bint is_ambiguous

    cdef long total_records = 0
    cdef long primary_mapped_records = 0
    cdef long below_mapq = 0
    cdef long multimapping_discarded = 0
    cdef long no_annotated_contig = 0
    cdef long no_feature = 0
    cdef long ambiguous = 0
    cdef long assigned_records = 0
    cdef double assigned_weight = 0.0

    bam = pysam.AlignmentFile(bam_path, "rb")

    try:
        for contig, contig_index in gene_index.items():
            tid = bam.get_tid(contig)

            if tid >= 0:
                index_by_tid[tid] = contig_index

        for aln in bam.fetch(until_eof=True):
            total_records += 1

            if aln.is_unmapped:
                continue

            if aln.is_secondary:
                continue

            if aln.is_supplementary:
                continue

            primary_mapped_records += 1

            if aln.mapping_quality < mapq_min:
                below_mapq += 1
                continue

            contig_index = index_by_tid.get(aln.reference_id)

            if contig_index is None:
                no_annotated_contig += 1
                continue

            weight = _read_weight_fast(
                aln,
                multi_mode,
                use_nh_tag,
            )

            if weight == 0.0:
                multimapping_discarded += 1
                continue

            result = _assign_alignment_to_gene_fast(
                aln,
                contig_index,
                strand_mode,
            )

            gene_id = result[0]
            is_ambiguous = result[1]

            if is_ambiguous:
                ambiguous += 1
                continue

            if gene_id is None:
                no_feature += 1
                continue

            counts[gene_id] = counts[gene_id] + weight
            assigned_records += 1
            assigned_weight += weight

    finally:
        bam.close()

    return counts, (
        total_records,
        primary_mapped_records,
        below_mapq,
        multimapping_discarded,
        no_annotated_contig,
        no_feature,
        ambiguous,
        assigned_records,
        assigned_weight,
    )
