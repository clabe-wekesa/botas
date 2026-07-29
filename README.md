# BOTAS

## Bacterial Operon-Aware Transcriptome Alignment System

BOTAS is an integrated Python framework for bacterial RNA-seq analysis that combines read alignment, gene quantification, and operon inference within a single command-line toolkit. Unlike general-purpose RNA-seq pipelines that require multiple external programs, BOTAS provides a unified workflow specifically designed for bacterial transcriptomes, with native support for circular genomes, paired-end sequencing, strand-specific expression analysis, and operon-aware transcriptomics.

The framework is designed for reproducible bacterial RNA-seq analysis while remaining modular enough for integration into automated bioinformatics workflows.

---

# Features

### Alignment

* Native alignment engine optimized for bacterial genomes
* Single-end and paired-end RNA-seq alignment
* Native support for circular chromosomes and plasmids
* Automatic handling of reads spanning the origin of replication
* Edit-distance alignment using Edlib
* Mapping quality estimation
* BAM output compatible with standard downstream tools

### Gene Quantification

* Read and fragment counting from coordinate-sorted BAM files
* Accurate paired-end fragment reconstruction
* Strand-specific counting
* Automatic handling of overlapping paired reads
* TPM and RPKM normalization
* Quantification from one or multiple BAM files
* Summary statistics for assigned, ambiguous, unmapped and unassigned fragments

### Operon Inference

* De novo operon prediction from RNA-seq alignments
* Integrates genomic organization and transcriptional evidence
* Uses intergenic distance, strand consistency, expression continuity and paired-end connectivity
* Consensus operon detection across multiple samples
* TSV and GFF output formats
* Configurable confidence thresholds

### General

* Pure Python implementation
* Minimal external dependencies
* Modular architecture for extension and development
* Command-line interface designed for reproducible analyses

---

# Why BOTAS?

Most RNA-seq analysis workflows were originally developed for eukaryotic transcriptomes and therefore emphasize splice-aware alignment, exon structure and transcript reconstruction. Bacterial transcriptomes differ fundamentally because genes are densely organized, transcription frequently occurs as polycistronic operons, chromosomes are commonly circular, and strand-specific transcription is often essential for accurate expression analysis.

BOTAS addresses these challenges through a workflow developed specifically for bacterial RNA-seq. Alignment, quantification and operon inference are performed within a single framework, eliminating the need to combine multiple independent tools while maintaining transparent and reproducible analyses.

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

# Quick Start

## Align paired-end reads

```bash
botas align \
    -r reference.fasta \
    -1 reads_R1.fastq \
    -2 reads_R2.fastq \
    -o aligned.bam
```

## Align single-end reads

```bash
botas align \
    -r reference.fasta \
    -U reads.fastq \
    -o aligned.bam
```

## Quantify gene expression

```bash
botas quantify \
    -b aligned.bam \
    -g genes.gff \
    -o gene_counts.tsv
```

BOTAS reports raw gene counts together with normalized expression estimates and fragment assignment statistics.

## Infer operons

```bash
botas getOperons \
    -b aligned.bam \
    -g genes.gff \
    -o operons.tsv
```

Operons may also be exported directly as GFF annotations for downstream genome browsers and comparative analyses.

---

# Project Structure

```text
botas/
├── cli/          Command-line interface
├── core/         Alignment engine
├── data/         Internal reference resources
├── io/           FASTA, FASTQ and BAM utilities
├── operons/      Operon inference
├── quantify/     Gene quantification
└── rrna/         rRNA filtering
```

---

# Requirements

* Python 3.10 or newer
* pysam
* edlib
* biopython

Optional:

* tqdm (progress reporting)

---

# Documentation

```bash
botas --help
botas align --help
botas quantify --help
botas getOperons --help
```

---

# Citation

If you use BOTAS in published research, please cite:

> Wekesa, C. S. **BOTAS: Bacterial Operon-Aware Transcriptome Alignment System.**

Citation details will be updated following publication.

---

# License

BOTAS is distributed under the MIT License.

---

# Author

**Clabe Simiyu Wekesa**

GitHub: https://github.com/clabe-wekesa
