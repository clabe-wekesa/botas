# BOTAS

## An Integrated Bacterial RNA-seq Analysis Framework with Circular-Aware Alignment and Operon Inference

**Clabe Simiyu Wekesa**, Kelvin Kiprotich, John Muoma, Axel Mithöfer

BOTAS (Bacterial Operon-Aware Transcriptome Alignment System) is an integrated framework for bacterial RNA-seq analysis that combines reference indexing, read alignment, gene quantification, and operon inference within a single command-line application.

Unlike conventional RNA-seq workflows that require multiple independent software packages, BOTAS provides an end-to-end workflow designed specifically for bacterial transcriptomics. The framework supports circular chromosomes and plasmids, single-end and paired-end sequencing, fragment-level gene quantification, and RNA-seq-guided operon inference while maintaining compatibility with standard genomic file formats.

BOTAS is implemented in Python and is designed for reproducible, modular, and high-throughput bacterial transcriptomic analyses.

---

# Features

## Native Reference Indexing

* Native minimizer-based reference indexing
* Reusable BOTAS index format (`*.botas.idx`)
* Configurable k-mer and minimizer window sizes
* Support for linear and circular bacterial genomes
* Circular overhang indexing for reads spanning the origin
* Selective circular treatment of chromosomes or plasmids

## RNA-seq Alignment

* Native bacterial read alignment engine
* Single-end and paired-end read alignment
* Support for circular genome boundary crossings
* Edit-distance alignment using Edlib
* Mapping-quality estimation
* Optional rRNA read filtering
* Optional coordinate sorting and BAM indexing
* Multiprocessing support

## Gene Quantification

* Gene-level paired-end fragment counting
* Reconstruction of paired-end fragments from coordinate-sorted BAM files
* Multi-sample quantification
* Configurable minimum feature overlap
* Optional largest-overlap assignment
* Assignment statistics for mapped, ambiguous, unmapped, and unassigned fragments
* Gene counts compatible with featureCounts under equivalent settings

## Operon Inference

* RNA-seq-guided operon prediction
* Integration of genomic adjacency and transcriptional evidence
* Strand-consistency analysis
* Intergenic-distance evaluation
* Gene-coverage comparison
* Consensus operon inference across multiple BAM files
* TSV output
* Optional GFF output

## General

* Python implementation
* Standard FASTA, FASTQ, BAM, GFF3, and TSV file support
* Reproducible command-line workflow
* Modular architecture
* Suitable for integration into automated pipelines

---

# Installation

## Install from PyPI

```bash
pip install botas-rnaseq
```

## Install the latest development version

```bash
pip install git+https://github.com/clabe-wekesa/botas.git
```

## Install from source

```bash
git clone https://github.com/clabe-wekesa/botas.git
cd botas
pip install .
```

## Development installation

```bash
pip install -e ".[dev]"
```

---

# Workflow

A typical BOTAS analysis consists of four steps.

```text
Reference FASTA
       │
       ▼
botas index
       │
       ▼
BOTAS index (.botas.idx)
       │
       ▼
botas align
       │
       ▼
Coordinate-sorted BAM
       │
       ├──────────────────┐
       ▼                  ▼
botas quantify     botas getOperons
       │                  │
       ▼                  ▼
 Gene counts       Operon predictions
```

---

# Quick Start

## 1. Build a reference index

```bash
botas index \
    --ref reference.fasta \
    --out reference.botas.idx
```

For a circular bacterial genome:

```bash
botas index \
    --ref reference.fasta \
    --circular \
    --out reference.botas.idx
```

Selected contigs, such as plasmids, can be treated as circular:

```bash
botas index \
    --ref reference.fasta \
    --circular-contigs plasmid1,plasmid2 \
    --out reference.botas.idx
```

When `--out` is omitted, BOTAS generates an output name ending in `.botas.idx`.

---

## 2. Align sequencing reads

### Paired-end alignment

```bash
botas align \
    --index reference.botas.idx \
    --fq1 reads_R1.fastq.gz \
    --fq2 reads_R2.fastq.gz \
    --pool \
    --threads 8 \
    --sort-bam \
    --out sample.bam
```

### Single-end alignment

```bash
botas align \
    --index reference.botas.idx \
    --fq reads.fastq.gz \
    --threads 4 \
    --sort-bam \
    --out sample.bam
```

BOTAS can also align directly from a reference FASTA:

```bash
botas align \
    --ref reference.fasta \
    --fq1 reads_R1.fastq.gz \
    --fq2 reads_R2.fastq.gz \
    --circular \
    --pool \
    --threads 8 \
    --sort-bam \
    --out sample.bam
```

The `--sort-bam` option creates a coordinate-sorted BAM file and its corresponding `.bai` index.

### Optional rRNA filtering

```bash
botas align \
    --index reference.botas.idx \
    --fq1 reads_R1.fastq.gz \
    --fq2 reads_R2.fastq.gz \
    --filter-rrna \
    --pool \
    --threads 8 \
    --sort-bam \
    --out sample.bam
```

BOTAS uses its bundled rRNA database unless a custom FASTA file is supplied with `--rrna-db`.

### Logging level

```bash
botas align ... --log INFO
botas align ... --log DEBUG
botas align ... --log WARN
botas align ... --log ERROR
```

The logging level must be provided after `--log`.

---

## 3. Quantify gene expression

```bash
botas quantify \
    --bam sample.sorted.bam \
    --gff annotation.gff \
    --feature-type gene \
    --id-attribute locus_tag \
    --out sample.gene_counts.tsv
```

The BAM input must contain coordinate-sorted paired-end alignments.

When `--id-attribute` is omitted, BOTAS uses `locus_tag` for gene features.

### Quantify multiple BAM files

```bash
botas quantify \
    --bam sample1.sorted.bam sample2.sorted.bam sample3.sorted.bam \
    --gff annotation.gff \
    --feature-type gene \
    --id-attribute locus_tag \
    --out gene_expression_matrix.tsv
```

### Largest-overlap assignment

By default, fragments overlapping multiple features remain ambiguous. They may instead be assigned to the feature with the largest overlap:

```bash
botas quantify \
    --bam sample.sorted.bam \
    --gff annotation.gff \
    --largest-overlap
```

Equal largest overlaps remain ambiguous.

---

## 4. Infer operons

### Infer operons from one BAM file

```bash
botas getOperons \
    --bam sample.sorted.bam \
    --gff annotation.gff \
    --out sample.operons.tsv
```

### Infer consensus operons from multiple BAM files

```bash
botas getOperons \
    --bam sample1.sorted.bam sample2.sorted.bam sample3.sorted.bam \
    --gff annotation.gff \
    --consensus \
    --out consensus.operons.tsv
```

### Write operons as GFF features

```bash
botas getOperons \
    --bam sample.sorted.bam \
    --gff annotation.gff \
    --write-gff \
    --prefix sample
```

Operon inference can be controlled using options such as:

* `--max-igd`
* `--min-coverage`
* `--min-cov-ratio`
* `--min-support`
* `--min-score`

---

# Working Directory

BOTAS creates a working directory for each analysis. A custom location can be specified using:

```bash
botas -d sample.botas align ...
```

The directory contains analysis logs, temporary files, final results, and a run manifest.

A typical structure is:

```text
sample.botas/
├── logs/
├── tmp/
├── results/
│   ├── sample.bam
│   ├── sample.sorted.bam
│   ├── sample.sorted.bam.bai
│   ├── sample.gene_counts.tsv
│   └── sample.operons.tsv
└── manifest.json
```

The exact files depend on the command and options used.

---

# Default Output Names

When `--out` is omitted, BOTAS generates an output name from the input filename.

| Command                             | Default suffix                  |
| ----------------------------------- | ------------------------------- |
| `botas index`                       | `.botas.idx`                    |
| `botas align`                       | `.bam`                          |
| Single-sample gene quantification   | `.gene_counts.tsv`              |
| Single-sample operon quantification | `.operon_counts.tsv`            |
| Multi-sample gene quantification    | `.gene_expression_matrix.tsv`   |
| Multi-sample operon quantification  | `.operon_expression_matrix.tsv` |
| `botas getOperons`                  | `.operons.tsv`                  |

BOTAS also appends the appropriate suffix when an output prefix is supplied without the expected extension.

---

# Command Overview

```text
botas index         Build a BOTAS reference index
botas align         Align RNA-seq reads
botas quantify      Quantify gene or operon expression
botas getOperons    Infer bacterial operons
```

Detailed help is available for every command:

```bash
botas --help
botas index --help
botas align --help
botas quantify --help
botas getOperons --help
```

The installed version can be displayed with:

```bash
botas --version
```

---

# Supported Input Formats

| Analysis         | Input                                     |
| ---------------- | ----------------------------------------- |
| Indexing         | FASTA                                     |
| Alignment        | FASTQ or compressed FASTQ                 |
| Quantification   | Coordinate-sorted paired-end BAM and GFF3 |
| Operon inference | Coordinate-sorted BAM and GFF3            |

---

# Supported Output Formats

| Analysis         | Output                      |
| ---------------- | --------------------------- |
| Indexing         | BOTAS index (`*.botas.idx`) |
| Alignment        | BAM and optional BAI        |
| Quantification   | TSV                         |
| Operon inference | TSV and optional GFF        |

---

# Requirements

* Python 3.10 or later
* pysam
* edlib
* Biopython

---

# Citation

When using BOTAS in published research, please cite:

> Wekesa CS, Kiprotich K, Muoma J, Mithöfer A.
> **BOTAS: An Integrated Bacterial RNA-seq Analysis Framework with Circular-Aware Alignment and Operon Inference.**
> Manuscript under review.

Citation information will be updated following publication.

---

# License

BOTAS is distributed under the MIT License.

---

# Author

**Clabe Simiyu Wekesa**

GitHub: https://github.com/clabe-wekesa
