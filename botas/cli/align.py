#!/usr/bin/env python3
# botas/cli/align.py

from __future__ import annotations

import argparse
import logging
import sys
import pysam
from pathlib import Path
from dataclasses import replace
from importlib.metadata import PackageNotFoundError, version
from botas.operons.cli import add_operon_args
from botas.core.align_core import align_read
from botas.core.parallel_context import PEContext
from botas.core.pe_pool import align_paired_pool
from botas.core.se_pool import align_single_pool
from botas.core.ref_index import KmerIndex
from botas.io.bam_writer import open_bam_writer, write_hit, write_pair, write_unmapped
from botas.io.fastq import read_fastq
from botas.io.reference_set import load_reference_set
from botas.operons.cli import run_get_operons
from botas.quantify.cli import add_quantify_args
from botas.io.workdir import init_workdir, update_manifest
from botas.operons.cli import GETOPERONS_EPILOG
from botas.quantify.cli import QUANTIFY_EPILOG
from botas.rrna.rrna import load_rrna_kmers, is_rrna_like
from botas.rrna.db import get_default_rrna_db
from botas.core.utils import normalize_read_name  # FIX: use fast name splitter
from botas.core.botas_index import build_botas_index, save_botas_index, default_index_path, BOTAS_INDEX_SUFFIX
from botas.core.botas_index import load_botas_index, BotasIntIndexAdapter
from botas.core.pe_pool import _normalize_circular_hit


ALIGN_EPILOG = """
Examples
--------
1) Paired-end alignment to a circular bacterial genome

  botas align --ref genome.fna --fq1 reads_R1.fq.gz --fq2 reads_R2.fq.gz --circular --threads 6 -o sample.bam


2) Single-end alignment

  botas align --ref genome.fna --fq reads.fq.gz --threads 4 -o sample.bam


3) Treat only selected contigs as circular (e.g. plasmids)

  botas align --ref genome.fna --fq1 reads_R1.fq.gz --fq2 reads_R2.fq.gz --circular-contigs plasmid1,plasmid2 -o sample.bam


Notes
-----
- Either --fq (single-end) OR (--fq1, --fq2) (paired-end) must be provided.
- Circular alignment enables correct handling of reads spanning the origin.
- Output BAM files are coordinate-sorted and suitable for downstream analysis.
"""


# ----------------------------
# Version / banner
# ----------------------------

def get_botas_version() -> str:
    try:
        return version("botas-rnaseq")
    except PackageNotFoundError:
        return "unknown"


def botas_banner() -> None:
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RESET = "\033[0m"
    BOLD = "\033[1m"

    banner = f"""
{CYAN}=================================================================={RESET}
{BOLD}{GREEN}                           BOTAS v{get_botas_version()}{RESET}
{YELLOW}              Bacterial Operon-aware Transcriptome Aligner Software{RESET}
{BOLD}{GREEN}      Developed by Dr. Clabe Wekesa at the MPI-CE (Jena, Germany){RESET}
{GREEN}               in Dr. Axel Mithöfer's lab{RESET}
{CYAN}=================================================================={RESET}
"""
    print(banner)


# ----------------------------
# Argparse helpers
# ----------------------------

class botasArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"[botas] ERROR: {message}\n")
        sys.stderr.write("For help: botas -h/--help\n")
        sys.exit(2)

    def print_usage(self, file=None):
        return  # suppress default argparse "usage:" spam


def _normalize_pair_name(name: str) -> str:
    """
    Return the read base name, stripping /1 or /2 suffix and any trailing
    fields after the first whitespace.

    FIX: was name.split()[0] — allocates a full list for every read.
    Now uses str.partition which stops at the first separator.
    """
    n = normalize_read_name(name)  # fast partition-based split
    if n.endswith("/1") or n.endswith("/2"):
        n = n[:-2]
    return n


def read_fastq_pairs(fq1: str, fq2: str):
    it1 = iter(read_fastq(fq1))
    it2 = iter(read_fastq(fq2))

    while True:
        r1 = next(it1, None)
        r2 = next(it2, None)

        if r1 is None and r2 is None:
            return
        if r1 is None or r2 is None:
            raise ValueError("FASTQ files have different number of records")

        n1 = _normalize_pair_name(r1.name)
        n2 = _normalize_pair_name(r2.name)
        if n1 != n2:
            raise ValueError(f"Mate name mismatch: {r1.name} vs {r2.name}")

        yield replace(r1, name=n1), replace(r2, name=n1)


class CleanHelpFormatter(argparse.RawTextHelpFormatter):
    """Custom formatter that shows -t/--threads without the THREADS metavar."""
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)

        parts = '/'.join(action.option_strings)

        if action.metavar is not None and action.metavar != '':
            parts += ' ' + self._format_args(action, action.metavar)
        elif action.nargs is not None and action.nargs != 0:
            default_metavar = self._get_default_metavar_for_optional(action)
            if default_metavar:
                parts += ' ' + self._format_args(action, default_metavar)

        return parts

def add_index_args(p: argparse.ArgumentParser) -> None:
    io = p.add_argument_group("Input / Output")
    io.add_argument(
        "-r", "--ref",
        required=True,
        metavar="FASTA",
        help="Reference FASTA file to index.",
    )
    io.add_argument(
        "-o", "--out",
        default=None,
        metavar="INDEX",
        help=(
            "Output BOTAS index file. "
            f"If omitted, BOTAS writes <reference>{BOTAS_INDEX_SUFFIX}."
        ),
    )

    idx = p.add_argument_group("Index parameters")
    idx.add_argument(
        "-k", "--kmer",
        type=int,
        default=15,
        metavar="INT",
        help="K-mer size used for BOTAS seed indexing. Default: 15.",
    )
    idx.add_argument(
        "-w", "--window",
        type=int,
        default=10,
        metavar="INT",
        help="Minimizer window size. Default: 10.",
    )
    idx.add_argument(
        "--circular",
        action="store_true",
        help="Treat all contigs as circular while building the index.",
    )
    idx.add_argument(
        "--circular-contigs",
        default=None,
        metavar="LIST",
        help="Comma-separated list of contig names to treat as circular.",
    )

    idx.add_argument(
        "--circular-overhang-percent",
        type=int,
        default=5,
        metavar="0-50",
        help=(
            "Percent of each circular contig end to append as overhang. "
            "0=no circular padding. "
            "5=recommended default. "
            "50=full circular doubling equivalent. "
            "Default: 5."
        ),
    )

    p.set_defaults(func=run_index)


def run_index(args: argparse.Namespace) -> int:
    circular_contigs = None
    if args.circular_contigs:
        circular_contigs = [
            x.strip()
            for x in args.circular_contigs.split(",")
            if x.strip()
        ]

    out = args.out or default_index_path(args.ref)

    idx = build_botas_index(
        args.ref,
        k=args.kmer,
        w=args.window,
        circular=args.circular,
        circular_contigs=circular_contigs,
        circular_overhang_percent=args.circular_overhang_percent,
    )

    save_botas_index(idx, out)

    total_contigs = len(idx.contigs)
    total_bases = sum(c.orig_len for c in idx.contigs)
    total_seeds = sum(len(c.index) for c in idx.contigs)
    total_positions = sum(sum(len(v) for v in c.index.values()) for c in idx.contigs)

    print(f"[botas index] reference      : {args.ref}")
    print(f"[botas index] output         : {out}")
    print(f"[botas index] contigs        : {total_contigs}")
    print(f"[botas index] bases          : {total_bases}")
    print(f"[botas index] k              : {args.kmer}")
    print(f"[botas index] w              : {args.window}")

    if args.circular:
        print(f"[botas index] circ overhang %: {args.circular_overhang_percent}")
    else:
        print("[botas index] circ overhang %: disabled (linear reference)")

    print(f"[botas index] unique seeds   : {total_seeds}")
    print(f"[botas index] seed positions : {total_positions}")
    print("[botas index] done")

    return 0

def build_root_parser() -> argparse.ArgumentParser:

    p = botasArgumentParser(
        prog="botas",
        description=(
            "BOTAS — a seed-and-extend aligner optimized for bacterial genomes, "
            "with native support for circular chromosomes and plasmids."
        ),
        formatter_class=CleanHelpFormatter,
        allow_abbrev=False,
    )

    p.add_argument(
        "--dir", "-d",
        dest="workdir",
        default=None,
        metavar="DIR",
        help=(
            "BOTAS working directory. "
            "If omitted, a directory is auto-created "
            "(e.g. botas_outdir/ or botas_run_YYYYmmdd_HHMMSS)."
        ),
    )

    p.add_argument("--version", action="version", version=f"BOTAS v{get_botas_version()}", help="Print version and exit.")

    sub = p.add_subparsers(dest="command", metavar="<command>", required=True)

    # ---------------- index ----------------
    index_p = sub.add_parser(
        "index",
        help="Build a BOTAS-native reference index",
        description=(
            "Build a BOTAS-native seed index for bacterial genomes. "
            "The output file uses the suffix .botas.idx."
        ),
        formatter_class=CleanHelpFormatter,
    )
    add_index_args(index_p)

    # ---------------- align ----------------
    align_p = sub.add_parser(
        "align",
        help="Align reads to a reference genome",
        description=(
            "Seed-and-extend alignment of RNA-seq reads to bacterial genomes, "
            "with native support for circular chromosomes and plasmids."
        ),
        formatter_class=CleanHelpFormatter,
        epilog=ALIGN_EPILOG,
    )
    add_align_args(align_p)

    # ---------------- getOperons ----------------
    operon_p = sub.add_parser(
        "getOperons",
        help="Infer bacterial operons from BAM alignments",
        description=(
            "Infer bacterial operons from RNA-seq alignments using "
            "genomic adjacency, strand consistency, and coverage coherence."
        ),
        formatter_class=CleanHelpFormatter,
        epilog=GETOPERONS_EPILOG,
    )
    add_operon_args(operon_p)

    # ---------------- quantify ----------------
    quantify_p = sub.add_parser(
        "quantify",
        help="Quantify gene or operon expression from BAM alignments",
        formatter_class=CleanHelpFormatter,
        epilog=QUANTIFY_EPILOG,
    )
    add_quantify_args(quantify_p)

    return p


# ----------------------------
# Parsers
# ----------------------------

def add_align_args(al: argparse.ArgumentParser) -> None:
    # =========================================================
    # Core input / output
    # =========================================================
    io0 = al.add_argument_group("Input / Output")
    io0.add_argument("-r", "--ref", required=False, metavar="FASTA", help="Reference FASTA file. Required unless --index is provided.")
    io0.add_argument("--index", default=None, metavar="BOTAS_IDX", help="BOTAS-native reference index file produced by 'botas index'.")
    io0.add_argument("-o", "--out", required=True, metavar="BAM", help="Output BAM file (overwritten if it exists).")

    # =========================================================
    # Read input
    # =========================================================
    reads = al.add_argument_group("Read input")
    reads.add_argument("--fq", metavar="FASTQ", help="Single-end FASTQ (.fq/.fastq, optionally .gz).")
    reads.add_argument("-1", "--fq1", metavar="FASTQ", help="Paired-end FASTQ read 1 (R1).")
    reads.add_argument("-2", "--fq2", metavar="FASTQ", help="Paired-end FASTQ read 2 (R2).")

    # =========================================================
    # Reference & seeding
    # =========================================================
    ref = al.add_argument_group("Reference and seeding")
    ref.add_argument("-c", "--circular", action="store_true", help="Treat all reference contigs as circular.")
    ref.add_argument("--circular-contigs", metavar="NAMES,NAMES", type=lambda s: s.split(","), default=None,
        help=(
            "Comma-separated contig names to treat as circular "
            "(overrides --circular for selected contigs)."))
    ref.add_argument(
        "--circular-overhang-percent",
        type=int,
        default=5,
        metavar="0-50",
        help=(
            "Percent of each circular contig end to append as overhang. "
            "0=no circular padding. "
            "5=recommended default. "
            "50=full circular doubling equivalent. "
            "Default: 5."
        ),
    )
    ref.add_argument("-k", "--kmer", type=int, default=15, help="Seed k-mer length (default: 15).")
    ref.add_argument("--step", type=int, default=5, help="Seed step size (default: 5).")
    ref.add_argument("--pad", type=int, default=250, help="Padding around alignment window in bp (default: 250).")
    ref.add_argument("--max-windows", type=int, default=50, help="Maximum candidate windows per read (default: 50).")
    ref.add_argument("--min-seed-hits", type=int, default=1,
        help="Minimum seed hits required to attempt DP alignment (default: 1).")

    # =========================================================
    # Paired-end geometry
    # =========================================================
    pe = al.add_argument_group("Paired-end geometry")
    pe.add_argument("--max-insert", type=int, default=1200, metavar="BP", help="Maximum allowed insert size in bp (default: 1200).")
    pe.add_argument("--expected-insert", type=int, default=300, metavar="BP", help="Expected insert size for mate rescue (default: 300).")
    pe.add_argument("--rescue-pad", type=int, default=200, metavar="BP", help="Search radius for mate rescue in bp (default: 200).")
    pe.add_argument("--no-rescue", action="store_true", help="Disable mate rescue (faster, less sensitive).")

    # =========================================================
    # Accuracy and debugging
    # =========================================================
    acc = al.add_argument_group("Accuracy and debugging")
    acc.add_argument("--sensitive", action="store_true", help="Enable more sensitive alignment mode (slower).")
    acc.add_argument("--debug-pairs", type=int, default=0, help="Debug the first N read pairs (paired-end only).")

    # =========================================================
    # Performance
    # =========================================================
    perf = al.add_argument_group("Performance")
    perf.add_argument("-t", "--threads", type=int, default=1, help="Number of worker processes (default: 1).")
    perf.add_argument("--pool", action="store_true", help="Use multiprocessing pool for alignment.")
    perf.add_argument("--chunk-size", type=int, default=2000, help="Chunk size per task when using --pool (default: 2000).")

    # =========================================================
    # BAM post-processing
    # =========================================================
    post = al.add_argument_group("BAM post-processing")
    post.add_argument("--sort-bam", action="store_true", help="Sort output BAM by coordinate and create index (.bai).")
    post.add_argument("--sort-threads", type=int, default=None, help="Threads for BAM sorting (default: same as --threads).")

    # =========================================================
    # rRNA filtering
    # =========================================================
    rrna = al.add_argument_group("rRNA filtering")
    rrna.add_argument("--filter-rrna", action="store_true",
                      help="Filter rRNA-like reads prior to alignment (k-mer based, fast).")
    rrna.add_argument("--rrna-db", metavar="FASTA",
                      help=("rRNA reference FASTA (.fa or .fa.gz). "
                            "If not provided, a curated bundled database is used."))
    rrna.add_argument("--rrna-k", type=int, default=18, metavar="K",
                      help="k-mer size for rRNA detection (default: 17).")
    rrna.add_argument("--rrna-min-hits", type=int, default=50, metavar="N",
                      help="Minimum shared rRNA k-mers to classify as rRNA (default: 3).")

    # =========================================================
    # Logging
    # =========================================================
    al.add_argument("--log", default="INFO", choices=["DEBUG", "INFO", "WARN", "ERROR"], help="Log verbosity level.")


# ----------------------------
# Command implementations
# ----------------------------

def run_align(args) -> int:
    # 1) Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(args._logfile, mode="a"), logging.StreamHandler(sys.stderr)],
    )
    log = logging.getLogger("botas.cli.align")

    # -------------------------------------------------
    # Optional rRNA filtering setup
    # -------------------------------------------------
    rrna_kmers = None

    if args.filter_rrna:
        rrna_db = args.rrna_db
        if rrna_db is None:
            rrna_db = str(get_default_rrna_db())
            log.info(f"[rrna] using bundled rRNA DB: {rrna_db}")
        else:
            log.info(f"[rrna] using user rRNA DB: {rrna_db}")

        rrna_kmers = load_rrna_kmers(rrna_db, k=args.rrna_k)

    rrna_examined = 0
    rrna_filtered = 0

    # 2) Validate inputs (SE xor PE)
    se = bool(args.fq)
    pe = bool(args.fq1 or args.fq2)
    if se and pe:
        log.error("Choose either --fq (SE) OR (--fq1 --fq2) (PE), not both.")
        return 2
    if not se and not (args.fq1 and args.fq2):
        log.error("Provide --fq (SE) or both --fq1 and --fq2 (PE).")
        return 2
    if (args.fq1 and not args.fq2) or (args.fq2 and not args.fq1):
        log.error("Paired-end requires both --fq1 and --fq2.")
        return 2

    # 3) Banner
    if log.isEnabledFor(logging.INFO):
        botas_banner()

    # 4) Load reference/index
    using_botas_index = args.index is not None

    if using_botas_index:
        idx_obj = load_botas_index(args.index)

        class _IndexedRefSet:
            def __init__(self, contigs):
                self._contigs = contigs

            def contigs(self):
                return self._contigs

            def __len__(self):
                return len(self._contigs)

        for c in idx_obj.contigs:
            c.index = BotasIntIndexAdapter(
                seq=c.seq,
                index=c.index,
                k=c.k,
                w=c.w,
                circular=c.circular,
                label=c.name,
            )

        refset = _IndexedRefSet(idx_obj.contigs)
        args.ref = idx_obj.ref_fa

        log.info("Loaded BOTAS index: %s", args.index)

    else:
        if not args.ref:
            log.error("Provide either --ref FASTA or --index BOTAS_IDX.")
            return 2

        refset = load_reference_set(
            args.ref,
            circular=args.circular,
            circular_contigs=getattr(args, "circular_contigs", None),
            circular_overhang_percent=args.circular_overhang_percent,
        )

    log.info("Loaded reference: %d contigs", len(refset))
    for c in refset.contigs():
        log.info("  - %s (len=%d) circular=%s", c.name, c.length, c.circular)

    # 5) Build per-contig indexes only when needed.
    # If --index is used, indexes are already loaded from disk.
    if not using_botas_index:
        if not (args.fq1 and args.fq2 and args.pool):
            for c in refset.contigs():
                c.index = KmerIndex(c.seq, k=args.kmer, circular=False, label=c.name)

    # 6) Open BAM writer (multi-contig header)
    bam_ref_lengths = []
    for c in refset.contigs():
        if c.circular:
            bam_ref_lengths.append(getattr(c, "orig_len", None) or (c.length // 2))
        else:
            bam_ref_lengths.append(c.length)

    bw = open_bam_writer(
        args.out,
        ref_names=[c.name for c in refset.contigs()],
        ref_lengths=bam_ref_lengths,
    )

    mapped = 0
    total = 0

    try:
        # =========================================================
        # PAIRED-END MODE
        # =========================================================
        if args.fq1 and args.fq2:
            log.info("Running paired-end mode: fq1=%s fq2=%s", args.fq1, args.fq2)

            ctx = PEContext(
                ref_fa=args.ref,
                circular_all=bool(args.circular),
                circular_overhang_percent=args.circular_overhang_percent,
                k=args.kmer,
                step=args.step,
                pad=args.pad,
                max_windows=args.max_windows,
                min_seed_hits=args.min_seed_hits,
                max_insert=args.max_insert,
                expected_insert=args.expected_insert,
                rescue_pad=args.rescue_pad,
                do_rescue=not args.no_rescue,
                sensitive=args.sensitive,
                debug_pairs=args.debug_pairs,
                prebuilt_contigs=list(refset.contigs()) if using_botas_index else None,
            )

            def _pe_iter_filtered(fq1, fq2):
                nonlocal rrna_examined, rrna_filtered

                for idx, (r1, r2) in enumerate(read_fastq_pairs(fq1, fq2)):
                    rrna_examined += 1

                    if rrna_kmers is not None:
                        if (
                            is_rrna_like(r1.seq, rrna_kmers, args.rrna_k, args.rrna_min_hits)
                            or
                            is_rrna_like(r2.seq, rrna_kmers, args.rrna_k, args.rrna_min_hits)
                        ):
                            rrna_filtered += 1
                            continue

                    yield idx, r1, r2

            pe_iter = _pe_iter_filtered(args.fq1, args.fq2)

            for idx, qname, r1s, r1q, r2s, r2q, ph in align_paired_pool(
                pairs_iter=pe_iter,
                ctx=ctx,
                threads=args.threads,
                chunk_size=args.chunk_size,
            ):
                total += 2
                if ph.hit1:
                    mapped += 1
                if ph.hit2:
                    mapped += 1

                write_pair(
                    bw,
                    qname=qname,
                    r1_seq=r1s,
                    r1_qual=r1q,
                    r2_seq=r2s,
                    r2_qual=r2q,
                    hit1=ph.hit1,
                    hit2=ph.hit2,
                    proper_pair=ph.proper_pair,
                    insert_size=ph.insert_size,
                )

                if (idx + 1) % 10000 == 0:
                    log.info("Processed %d pairs...", idx + 1)

            log.info(
                "Done (PE). total_reads=%d mapped_reads=%d (%.2f%%)",
                total, mapped, 100.0 * mapped / max(1, total),
            )

        # =========================================================
        # SINGLE-END MODE
        # =========================================================
        else:
            log.info("Running single-end mode: fq=%s", args.fq)

            def _se_iter_filtered(fq):
                nonlocal rrna_examined, rrna_filtered

                for idx, r in enumerate(read_fastq(fq)):
                    rrna_examined += 1

                    if rrna_kmers is not None:
                        if is_rrna_like(
                            r.seq,
                            rrna_kmers,
                            k=args.rrna_k,
                            min_hits=args.rrna_min_hits,
                        ):
                            rrna_filtered += 1
                            continue

                    yield idx, r

            reads_iter = _se_iter_filtered(args.fq)

            # pooled SE
            if args.pool and args.threads > 1:
                se_iter = align_single_pool(
                    reads_iter=reads_iter,
                    ref_fa=args.ref,
                    threads=args.threads,
                    chunk_size=args.chunk_size,
                    k=args.kmer,
                    step=args.step,
                    pad=args.pad,
                    max_windows=args.max_windows,
                    min_seed_hits=args.min_seed_hits,
                    prebuilt_contigs=list(refset.contigs()) if using_botas_index else None,
                )
            # serial SE
            else:
                def _serial_se():
                    for idx, r in reads_iter:
                        best_hit = None
                        second_hit = None

                        for c in refset.contigs():
                            hit = align_read(
                                read_seq=r.seq,
                                rname=c.name,
                                ref_seq=c.seq,
                                index=c.index,
                                circular=False,
                                k=args.kmer,
                                step=args.step,
                                pad=args.pad,
                                max_windows=args.max_windows,
                                min_seed_hits=args.min_seed_hits,
                            )
                            if hit is None:
                                continue

                            if best_hit is None or hit.ascore > best_hit.ascore:
                                second_hit = best_hit
                                best_hit = hit
                            elif second_hit is None or hit.ascore > second_hit.ascore:
                                second_hit = hit

                        yield idx, r.name, r.seq, r.qual, best_hit, second_hit

                se_iter = _serial_se()

            for idx, qname, seq, qual, best_hit, second_hit in se_iter:

                total += 1

                if best_hit is None:
                    write_unmapped(bw, qname=qname, seq=seq, qual=qual)
                else:
                    mapped += 1
                    if second_hit:
                        delta = best_hit.ascore - second_hit.ascore
                        mapq = max(0, min(60, 10 + delta))
                    else:
                        mapq = 60

                    best_hit = replace(best_hit, mapq=mapq)

                    # Convert padded circular coordinates back to original genome coordinates.
                    hit_contig = next(
                        cc for cc in refset.contigs()
                        if cc.name == best_hit.rname
                    )

                    best_hit = _normalize_circular_hit(best_hit, hit_contig)

                    write_hit(
                        bw,
                        qname=qname,
                        seq=seq,
                        qual=qual,
                        hit_rname=best_hit.rname,
                        pos0=best_hit.pos0,
                        cigar=best_hit.cigar,
                        mapq=best_hit.mapq,
                        strand=best_hit.strand,
                        nm=-best_hit.ascore,
                        ascore=best_hit.ascore,
                    )

                if (idx + 1) % 10000 == 0:
                    log.info("Processed %d reads...", idx + 1)

            log.info(
                "Done (SE). total=%d mapped=%d (%.2f%%)",
                total, mapped, 100.0 * mapped / max(1, total),
            )

    finally:
        bw.close()

    # -------------------------------------------------
    # Optional BAM sorting and indexing
    # -------------------------------------------------
    if getattr(args, "sort_bam", False):
        sort_threads = args.sort_threads or args.threads
        out_bam = Path(args.out)
        if out_bam.name.endswith(".sorted.bam"):
            sorted_bam = out_bam
        else:
            sorted_bam = out_bam.with_suffix(".sorted.bam")
        log.info("Sorting BAM (%d threads): %s → %s", sort_threads, out_bam.name, sorted_bam.name)
        pysam.sort("-@", str(sort_threads), "-o", str(sorted_bam), str(out_bam))
        log.info("Indexing BAM: %s", sorted_bam.name)
        pysam.index(str(sorted_bam))
        log.info("BAM sorting and indexing complete")

    # -------------------------------------------------
    # Report rRNA filtering statistics
    # -------------------------------------------------
    if rrna_kmers is not None:
        retained = rrna_examined - rrna_filtered
        pct = (rrna_filtered / rrna_examined * 100.0) if rrna_examined else 0.0
        unit = "pairs" if (args.fq1 and args.fq2) else "reads"
        log.info(
            "[rRNA] examined_%s=%d filtered_%s=%d (%.1f%%) retained_%s=%d",
            unit, rrna_examined,
            unit, rrna_filtered, pct,
            unit, retained,
        )

    return 0


# ----------------------------
# Main
# ----------------------------

def main(argv: list[str] | None = None) -> int:
    if argv is None and len(sys.argv) == 1:
        sys.stderr.write("[botas] ERROR: no command-line arguments provided\n")
        sys.stderr.write("For help: botas -h/--help\n")
        return 2

    parser = build_root_parser()
    args = parser.parse_args(argv)

    # -------------------------------------------------
    # Initialize BOTAS working directory
    # -------------------------------------------------
    wd = init_workdir(args.workdir, out=getattr(args, "out", None))

    # -------------------------------------------------
    # Attach workdir paths
    # -------------------------------------------------
    args._workdir = str(wd.root)
    args._tmpdir = str(wd.tmp)
    args._resultsdir = str(wd.results)
    args._logfile = str(wd.logfile)

    # -------------------------------------------------
    # Resolve output paths INSIDE workdir
    # -------------------------------------------------
    if hasattr(args, "out") and args.out is not None:
        out_name = Path(args.out).name
        args.out = str(Path(args._resultsdir) / out_name)

    # -------------------------------------------------
    # Update manifest — single call
    # FIX: was called twice in a row with overlapping keys; merged into one.
    # -------------------------------------------------
    update_manifest(
        wd,
        {
            "command": " ".join(sys.argv),
            "cwd": str(Path.cwd()),
            "workdir": str(wd.root),
            "resultsdir": args._resultsdir,
        },
    )

    if args.command is None:
        sys.stderr.write("[botas] ERROR: no command specified\n")
        sys.stderr.write("For help: botas -h/--help\n")
        return 2

    try:
        if args.command == "index":
            return args.func(args)

        if args.command == "align":
            return run_align(args)

        if args.command == "getOperons":
            return run_get_operons(args)

        if args.command == "quantify":
            return args.func(args)

        sys.stderr.write(f"[botas] ERROR: unknown command {args.command}\n")
        return 2

    except KeyboardInterrupt:
        logging.getLogger("botas").warning("⛔ Interrupted by user. Exiting.")
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
