#!/usr/bin/env python3

import argparse
import csv
import re
import sys
import time
from collections import deque
from urllib.parse import urljoin, urlparse, urldefrag

import requests
from bs4 import BeautifulSoup


# ============================================================
# SETTINGS
# ============================================================

DEFAULT_HOST = "dbtbs.hgc.jp"

ALL_OUT = None
EXPERIMENTAL_OUT = None
PAIR_OUT = None
FAILED_OUT = None

USER_AGENT = (
    "Mozilla/5.0 "
    "(compatible; BOTAS-operon-benchmark/1.0; "
    "+https://github.com/clabe-wekesa/botas)"
)

REQUEST_DELAY = 0.25


# ============================================================
# EXPERIMENTAL EVIDENCE
# ============================================================

EXPERIMENTAL_PATTERNS = [
    r"\bnorthern\b",
    r"\bnorthern blot",
    r"\brt[\s-]?pcr\b",
    r"\breverse transcription\b",
    r"\bprimer extension\b",
    r"\bs1 nuclease\b",
    r"\bs1 mapping\b",
    r"\btranscriptional fusion",
    r"\blacz\b",
    r"\breporter fusion",
    r"\brna analysis\b",
    r"\btranscript analysis\b",
    r"\btranscript size\b",
    r"\btranscription analysis\b",
    r"\b5['′]?\s*race\b",
    r"\brnase protection\b",
    r"\bhybridization\b",
]

NONEXPERIMENTAL_ONLY_PATTERNS = [
    r"^\s*genome analysis\s*$",
    r"^\s*prediction\s*$",
    r"^\s*computational prediction\s*$",
]


# ============================================================
# COMMAND LINE
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Scrape DBTBS regulated operons and construct an "
            "experimentally supported adjacent-gene pair reference."
        )
    )

    parser.add_argument(
        "--start-url",
        action="append",
        required=True,
        help=(
            "DBTBS index/category page containing links to operons. "
            "May be supplied more than once."
        ),
    )

    parser.add_argument(
        "--crawl-depth",
        type=int,
        default=1,
        help=(
            "Number of internal index-page levels to follow from each "
            "start page (default: 1)."
        ),
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=REQUEST_DELAY,
        help="Delay between HTTP requests in seconds (default: 0.25).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for DBTBS tables and failed-page log.",
    )

    return parser.parse_args()


# ============================================================
# HTTP
# ============================================================

session = requests.Session()

session.headers.update({
    "User-Agent": USER_AGENT,
})


def get_html(url, delay):

    print(f"GET {url}")

    response = session.get(
        url,
        timeout=30,
    )

    response.raise_for_status()

    if delay:
        time.sleep(delay)

    return response.text


def get_soup(url, delay):

    return BeautifulSoup(
        get_html(url, delay),
        "html.parser"
    )


# ============================================================
# URL HELPERS
# ============================================================

def normalize_url(url):

    url, _fragment = urldefrag(url)

    return url


def is_internal_dbtbs(url):

    parsed = urlparse(url)

    return parsed.netloc.lower() == DEFAULT_HOST


def is_operon_url(url):

    path = urlparse(url).path

    return bool(
        re.search(
            r"/(?:ver\d+/)?COG/prom/[^/]+\.html$",
            path,
            flags=re.I,
        )
    )


def likely_index_page(url):

    path = urlparse(url).path.lower()

    if is_operon_url(url):
        return False

    if path.endswith((
        ".jpg",
        ".jpeg",
        ".png",
        ".gif",
        ".css",
        ".js",
        ".pdf",
        ".gz",
        ".zip",
        ".txt",
    )):
        return False

    return True


# ============================================================
# DISCOVER OPERON LINKS
# ============================================================

def discover_operon_pages(start_urls, crawl_depth, delay):

    queue = deque()

    for url in start_urls:
        queue.append((normalize_url(url), 0))

    visited = set()
    operon_urls = set()

    while queue:

        url, depth = queue.popleft()

        if url in visited:
            continue

        visited.add(url)

        try:
            soup = get_soup(url, delay)

        except Exception as exc:
            print(
                f"WARNING: could not read index page {url}: {exc}",
                file=sys.stderr,
            )
            continue

        links_found_here = 0

        for link in soup.find_all("a", href=True):

            target = normalize_url(
                urljoin(url, link["href"])
            )

            if not is_internal_dbtbs(target):
                continue

            if is_operon_url(target):

                operon_urls.add(target)
                links_found_here += 1
                continue

            if (
                depth < crawl_depth
                and likely_index_page(target)
                and target not in visited
            ):
                queue.append(
                    (target, depth + 1)
                )

        print(
            f"  operon links found on page: "
            f"{links_found_here}"
        )

    print()
    print(f"Index pages visited: {len(visited)}")
    print(
        f"Unique operon pages discovered: "
        f"{len(operon_urls)}"
    )

    return sorted(operon_urls)


# ============================================================
# TEXT HELPERS
# ============================================================

def clean_text(value):

    if value is None:
        return ""

    return " ".join(
        value.replace("\xa0", " ").split()
    )


def classify_evidence(evidence):

    evidence_clean = clean_text(evidence)
    evidence_lower = evidence_clean.lower()

    if not evidence_lower:
        return False

    for pattern in NONEXPERIMENTAL_ONLY_PATTERNS:

        if re.search(
            pattern,
            evidence_lower,
            flags=re.I,
        ):
            return False

    for pattern in EXPERIMENTAL_PATTERNS:

        if re.search(
            pattern,
            evidence_lower,
            flags=re.I,
        ):
            return True

    return False


# ============================================================
# EXTRACT LABELLED FIELD
# ============================================================

def extract_labelled_field(soup, label):

    wanted = label.lower().rstrip(":")

    # --------------------------------------------------------
    # Tables
    # --------------------------------------------------------

    for row in soup.find_all("tr"):

        cells = row.find_all(
            ["th", "td"],
            recursive=False
        )

        if len(cells) < 2:
            continue

        first = clean_text(
            cells[0].get_text(" ", strip=True)
        )

        first_norm = first.lower().rstrip(":")

        if first_norm == wanted:

            return clean_text(
                cells[1].get_text(
                    " ",
                    strip=True
                )
            )

    # --------------------------------------------------------
    # Generic text fallback
    # --------------------------------------------------------

    text = soup.get_text(
        "\n",
        strip=True
    )

    lines = [
        clean_text(x)
        for x in text.splitlines()
        if clean_text(x)
    ]

    for i, line in enumerate(lines):

        normalized = (
            line.lower()
            .strip()
            .rstrip(":")
        )

        if normalized == wanted:

            if i + 1 < len(lines):
                return lines[i + 1]

        prefix = wanted + ":"

        if line.lower().startswith(prefix):

            return clean_text(
                line[len(prefix):]
            )

    return ""


# ============================================================
# EXTRACT GENES
# ============================================================

def extract_genes(soup):

    # --------------------------------------------------------
    # Strategy 1:
    # find tables whose header contains "Gene" or "Genes"
    # --------------------------------------------------------

    for table in soup.find_all("table"):

        rows = table.find_all("tr")

        if not rows:
            continue

        header_cells = rows[0].find_all(
            ["th", "td"]
        )

        headers = [
            clean_text(
                x.get_text(" ", strip=True)
            ).lower()
            for x in header_cells
        ]

        gene_col = None

        for i, header in enumerate(headers):

            if header in {
                "gene",
                "genes",
                "gene name",
            }:
                gene_col = i
                break

        if gene_col is None:
            continue

        genes = []

        for row in rows[1:]:

            cells = row.find_all(
                ["td", "th"]
            )

            if len(cells) <= gene_col:
                continue

            gene = clean_text(
                cells[gene_col].get_text(
                    " ",
                    strip=True
                )
            )

            if not gene:
                continue

            if gene.lower() in {
                "gene",
                "genes",
            }:
                continue

            genes.append(gene)

        if len(genes) >= 1:
            return genes

    # --------------------------------------------------------
    # Strategy 2:
    # some DBTBS pages represent genes as links to gene pages.
    # --------------------------------------------------------

    genes = []

    for a in soup.find_all("a", href=True):

        href = a["href"]

        text = clean_text(
            a.get_text(" ", strip=True)
        )

        if not text:
            continue

        if (
            "/gene/" in href.lower()
            or "/prom/" in href.lower()
        ):

            if re.fullmatch(
                r"[A-Za-z][A-Za-z0-9_.-]*",
                text
            ):
                genes.append(text)

    # remove consecutive duplicates
    cleaned = []

    for gene in genes:

        if not cleaned or cleaned[-1] != gene:
            cleaned.append(gene)

    return cleaned


# ============================================================
# OPERON NAME
# ============================================================

def get_operon_name(url, soup):

    filename = (
        urlparse(url)
        .path
        .rstrip("/")
        .split("/")[-1]
    )

    filename = re.sub(
        r"\.html$",
        "",
        filename,
        flags=re.I,
    )

    # Try page heading first

    for tag_name in ["h1", "h2", "h3"]:

        tag = soup.find(tag_name)

        if not tag:
            continue

        heading = clean_text(
            tag.get_text(" ", strip=True)
        )

        match = re.search(
            r"(?:regulated\s+)?operon\s*:?\s*(.+)",
            heading,
            flags=re.I,
        )

        if match:

            candidate = clean_text(
                match.group(1)
            )

            if candidate:
                return candidate

    return filename


# ============================================================
# PARSE OPERON PAGE
# ============================================================

def parse_operon_page(url, delay):

    soup = get_soup(url, delay)

    operon = get_operon_name(
        url,
        soup
    )

    genes = extract_genes(
        soup
    )

    evidence = extract_labelled_field(
        soup,
        "Operon evidence"
    )

    reference = extract_labelled_field(
        soup,
        "Reference"
    )

    experimental = classify_evidence(
        evidence
    )

    return {
        "operon": operon,
        "genes": genes,
        "n_genes": len(genes),
        "evidence": evidence,
        "experimental": experimental,
        "reference": reference,
        "url": url,
    }


# ============================================================
# WRITE OUTPUT
# ============================================================

def write_operons(records, filename):

    fields = [
        "operon",
        "genes",
        "n_genes",
        "evidence",
        "experimental",
        "reference",
        "url",
    ]

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for record in records:

            writer.writerow({
                "operon":
                    record["operon"],

                "genes":
                    ",".join(
                        record["genes"]
                    ),

                "n_genes":
                    record["n_genes"],

                "evidence":
                    record["evidence"],

                "experimental":
                    str(
                        record["experimental"]
                    ).lower(),

                "reference":
                    record["reference"],

                "url":
                    record["url"],
            })


# ============================================================
# ADJACENT OPERON PAIRS
# ============================================================

def build_truth_pairs(records):

    pairs = {}

    for record in records:

        genes = record["genes"]

        if len(genes) < 2:
            continue

        for i in range(
            len(genes) - 1
        ):

            gene1 = genes[i]
            gene2 = genes[i + 1]

            key = (
                gene1,
                gene2,
            )

            if key not in pairs:

                pairs[key] = {
                    "gene1": gene1,
                    "gene2": gene2,
                    "operons": set(),
                    "evidence": set(),
                    "urls": set(),
                }

            pairs[key]["operons"].add(
                record["operon"]
            )

            if record["evidence"]:
                pairs[key]["evidence"].add(
                    record["evidence"]
                )

            pairs[key]["urls"].add(
                record["url"]
            )

    return pairs


def write_pairs(pairs, filename):

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow([
            "gene1",
            "gene2",
            "operon",
            "evidence",
            "source_url",
        ])

        for key in sorted(pairs):

            record = pairs[key]

            writer.writerow([
                record["gene1"],
                record["gene2"],

                ";".join(
                    sorted(
                        record["operons"]
                    )
                ),

                "; ".join(
                    sorted(
                        record["evidence"]
                    )
                ),

                ";".join(
                    sorted(
                        record["urls"]
                    )
                ),
            ])


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    from pathlib import Path
    global ALL_OUT, EXPERIMENTAL_OUT, PAIR_OUT, FAILED_OUT
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    ALL_OUT = str(outdir / "dbtbs_operons_all.tsv")
    EXPERIMENTAL_OUT = str(outdir / "dbtbs_operons_experimental.tsv")
    PAIR_OUT = str(outdir / "dbtbs_operon_truth_pairs.tsv")
    FAILED_OUT = str(outdir / "dbtbs_failed_pages.tsv")

    print(
        "========================================"
    )
    print(
        "DBTBS OPERON SCRAPER"
    )
    print(
        "========================================"
    )

    print("Starting pages:")

    for url in args.start_url:
        print(f"  {url}")

    print()

    # --------------------------------------------------------
    # Discover operon pages
    # --------------------------------------------------------

    operon_urls = discover_operon_pages(
        args.start_url,
        args.crawl_depth,
        args.delay,
    )

    if not operon_urls:

        print()
        print(
            "ERROR: No operon pages were discovered."
        )
        print()
        print(
            "Use --start-url with the DBTBS page "
            "that contains the gene/operon tables."
        )

        sys.exit(1)

    # --------------------------------------------------------
    # Parse operons
    # --------------------------------------------------------

    records = []
    failures = []

    print()
    print(
        "========================================"
    )
    print(
        "PARSING OPERONS"
    )
    print(
        "========================================"
    )

    for index, url in enumerate(
        operon_urls,
        start=1,
    ):

        print(
            f"[{index}/{len(operon_urls)}]"
        )

        try:

            record = parse_operon_page(
                url,
                args.delay,
            )

        except Exception as exc:

            print(
                f"WARNING: failed {url}: {exc}",
                file=sys.stderr,
            )

            failures.append(
                (url, str(exc))
            )

            continue

        # keep the page even if only one gene;
        # singletons are useful for auditing
        records.append(record)

        print(
            f"  {record['operon']}: "
            f"{record['n_genes']} genes | "
            f"experimental="
            f"{record['experimental']}"
        )

    # --------------------------------------------------------
    # Remove duplicate operon pages
    #
    # DBTBS can expose the same operon from several genes.
    # --------------------------------------------------------

    deduplicated = {}

    for record in records:

        key = (
            record["operon"],
            tuple(record["genes"]),
        )

        if key not in deduplicated:
            deduplicated[key] = record

    records = list(
        deduplicated.values()
    )

    # --------------------------------------------------------
    # Strict experimental set
    # --------------------------------------------------------

    experimental = [
        record
        for record in records
        if (
            record["experimental"]
            and record["n_genes"] >= 2
        )
    ]

    # --------------------------------------------------------
    # Write outputs
    # --------------------------------------------------------

    write_operons(
        records,
        ALL_OUT,
    )

    write_operons(
        experimental,
        EXPERIMENTAL_OUT,
    )

    pairs = build_truth_pairs(
        experimental
    )

    write_pairs(
        pairs,
        PAIR_OUT,
    )

    # --------------------------------------------------------
    # Failures
    # --------------------------------------------------------

    with open(
        FAILED_OUT,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow([
            "url",
            "error",
        ])

        writer.writerows(
            failures
        )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print(
        "========================================"
    )
    print(
        "DBTBS SCRAPE SUMMARY"
    )
    print(
        "========================================"
    )

    print(
        f"Unique operon pages:          "
        f"{len(records)}"
    )

    print(
        f"Multi-gene experimental:     "
        f"{len(experimental)}"
    )

    print(
        f"Experimental adjacent pairs: "
        f"{len(pairs)}"
    )

    print(
        f"Failed pages:                 "
        f"{len(failures)}"
    )

    print()
    print("Outputs:")

    print(
        f"  {ALL_OUT}"
    )

    print(
        f"  {EXPERIMENTAL_OUT}"
    )

    print(
        f"  {PAIR_OUT}"
    )

    print(
        f"  {FAILED_OUT}"
    )


if __name__ == "__main__":
    main()
