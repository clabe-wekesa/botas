# BOTAS

## An Integrated Bacterial RNA-seq Analysis Framework with Circular-Aware Alignment and Operon Inference

**Clabe Simiyu Wekesa**, Kelvin Kiprotich, John Muoma, Axel Mithöfer

BOTAS (Bacterial Operon-Aware Transcriptome Alignment System) is an integrated framework for bacterial RNA-seq analysis that combines native reference indexing, read alignment, gene quantification, and operon inference within a single command-line application.

Unlike conventional RNA-seq workflows that require multiple independent software packages, BOTAS provides an end-to-end bacterial transcriptomics workflow specifically designed for prokaryotic genomes. The framework natively supports circular chromosomes and plasmids, paired-end sequencing, strand-specific expression analysis, and operon-aware transcriptome analysis while maintaining compatibility with standard genomic file formats.

BOTAS is implemented entirely in Python and is designed for reproducible, modular and high-throughput bacterial transcriptomic analyses.

---

# Features

## Native Reference Indexing

- Native minimizer-based reference indexing
- Reusable BOTAS index format (`*.botas.idx`)
- Configurable k-mer and minimizer window sizes
- Support for linear and circular bacterial genomes
- Circular overhang indexing for origin-spanning reads

## RNA-seq Alignment

- Native bacterial read alignment engine
- Single-end and paired-end read alignment
- Automatic handling of circular genome boundary crossings
- Edit-distance alignment using Edlib
- Mapping quality estimation
- Coordinate-sorted BAM output
- Multi-threaded execution

## Gene Quantification

- Gene-level read and fragment counting
- Accurate paired-end fragment reconstruction
- Strand-specific quantification
- TPM and RPKM normalization
- Multi-sample quantification
- Assignment statistics for mapped, ambiguous, unmapped and unassigned fragments
- Gene counts compatible with featureCounts

## Operon Inference

- RNA-seq-guided operon prediction
- Integration of genomic organization and transcriptional evidence
- Strand consistency analysis
- Intergenic distance evaluation
- Expression continuity assessment
- Paired-end connectivity analysis
- Consensus operon prediction across multiple samples
- TSV and GFF outputs

## General

- Pure Python implementation
- Minimal external dependencies
- Reproducible command-line workflow
- Modular architecture
- Easily integrated into automated pipelines

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

```
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
       ├──────────────┐
       ▼              ▼
botas quantify   botas getOperons
       │              │
       ▼              ▼
 Gene counts      Operon predictions
```

---

# Quick Start

## 1. Build a reference index

```bash
botas index -r reference.fasta -o reference.botas.idx
```

For circular bacterial genomes:

```bash
botas index -r reference.fasta -o reference.botas.idx --circular
```

---

## 2. Align sequencing reads

### Paired-end

```bash
botas align -x reference.botas.idx -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz --sort-bam
```

### Single-end

```bash
botas align -x reference.botas.idx -U reads.fastq.gz --sort-bam
```

BOTAS automatically creates a project directory containing intermediate files, logs and final alignment results.

---

## 3. Quantify gene expression

```bash
botas quantify -b alignment.bam -g annotation.gff --feature-type gene --gff-gene-attribute locus_tag
```

BOTAS reports

- raw gene counts
- TPM
- RPKM
- fragment assignment statistics
- summary reports

---

## 4. Infer operons

```bash
botas getOperons -b alignment.bam -g annotation.gff
```

Predicted operons are exported as

- TSV tables
- GFF annotations

for downstream visualization and comparative genomics analyses.

---

# Output Structure

Each BOTAS analysis creates a dedicated working directory.

```
sample.botas/
│
├── logs/
├── temp/
├── results/
│   ├── alignment.bam
│   ├── alignment.bam.bai
│   ├── gene_counts.tsv
│   ├── summary.tsv
│   ├── operons.tsv
│   └── operons.gff
│
└── config.json
```

---

# Command Overview

```
botas index         Build a BOTAS reference index

botas align         Align RNA-seq reads

botas quantify      Quantify gene expression

botas getOperons    Predict bacterial operons
```

Detailed help is available for every command.

```bash
botas --help
botas index --help
botas align --help
botas quantify --help
botas getOperons --help
```

---

# Supported Input Formats

| Analysis | Input |
|----------|-------|
| Indexing | FASTA |
| Alignment | FASTQ, FASTQ.GZ |
| Quantification | BAM, GFF3 |
| Operon inference | BAM, GFF3 |

---

# Supported Output Formats

| Analysis | Output |
|----------|--------|
| Indexing | BOTAS index (`*.botas.idx`) |
| Alignment | BAM |
| Quantification | TSV |
| Operon inference | TSV, GFF |

---

# Requirements

- Python 3.10 or later
- pysam
- edlib
- Biopython

Optional

- tqdm

---

# Citation

If you use BOTAS in published research, please cite:

> Wekesa CS, Kiprotich K, Muoma J, Mithöfer A.
> **BOTAS: An Integrated Bacterial RNA-seq Analysis Framework with Circular-Aware Alignment and Operon Inference.**
> *(Manuscript under review.)*

Citation information will be updated following publication.

---

# License

BOTAS is distributed under the MIT License.

---

# Author

**Clabe Simiyu Wekesa**

GitHub: https://github.com/clabe-wekesa
