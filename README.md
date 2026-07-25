# BOTAS

## Bacterial Operon-Aware Transcriptome Alignment System

BOTAS is a bacterial RNA-seq analysis framework that integrates read alignment, gene quantification, operon inference, and operon-level quantification within a unified Python package. It is specifically designed for bacterial transcriptomes, providing native support for circular genomes, operon-aware analyses, and strand-specific expression profiling.

Unlike general-purpose RNA-seq aligners developed primarily for eukaryotic transcriptomes, BOTAS addresses the unique characteristics of bacterial genomes, including dense gene organization, polycistronic transcription, operon architecture, and circular chromosomes.

---

## Key Features

- Seed-and-extend alignment engine optimized for bacterial genomes
- Native support for circular chromosomes and plasmids
- Paired-end and single-end read alignment
- Strand-aware gene expression quantification
- Gene-level and operon-level expression quantification
- Operon inference from aligned reads
- Operon-aware transcriptome analysis
- Optional rRNA filtering
- Parallel execution for scalable performance
- Modular and extensible Python architecture

---

## Why BOTAS?

Most RNA-seq aligners were designed for eukaryotic transcriptomes and assume splicing, large introns, and linear chromosomes. These assumptions are inappropriate for bacterial transcriptomes, which are characterized by:

- Densely packed genes
- Polycistronic transcripts organized into operons
- Circular chromosomes and plasmids
- Strong dependence on strand-specific transcription

BOTAS addresses these challenges through native circular-genome support, operon-aware analysis, strand-specific quantification, circular insert-size validation, and integrated gene and operon quantification within a single framework.

---

## Installation

### Install from PyPI

```bash
pip install botas
```

### Install the latest development version

```bash
pip install git+https://github.com/clabe-wekesa/botas.git
```

### Install from source

```bash
git clone https://github.com/clabe-wekesa/botas.git
cd botas
pip install .
```

### Development installation

```bash
pip install -e ".[dev]"
```

---

## Quick Start

### Paired-end alignment

```bash
botas align \
    -r reference.fasta \
    -1 reads_R1.fastq \
    -2 reads_R2.fastq \
    -o aligned.bam
```

### Single-end alignment

```bash
botas align \
    -r reference.fasta \
    -U reads.fastq \
    -o aligned.bam
```

### Gene quantification

```bash
botas quant \
    -b aligned.bam \
    -g genes.gff \
    -o gene_counts.tsv
```

Supported features include:

- Strand-specific counting
- MAPQ filtering
- Multi-mapper handling
- TPM and RPKM calculation
- Multi-BAM count matrix generation

### Operon inference

```bash
botas getOperons \
    -b aligned.bam \
    -g genes.gff \
    -o operons.tsv
```

Operon inference integrates:

- Intergenic distance
- Strand consistency
- Coverage similarity
- Paired-end support
- Consensus-based merging

---

## Architecture

```text
botas/
├── cli/          # Command-line interface
├── core/         # Alignment engine
├── data/         # Reference resources
├── io/           # FASTQ, BAM and reference handling
├── operons/      # Operon inference
├── quantify/     # Gene and operon quantification
└── rrna/         # rRNA detection and filtering
```

The alignment engine implements:

- K-mer indexing
- Seed clustering
- Edit-distance extension using edlib
- CIGAR reconstruction
- Mapping quality estimation
- Circular coordinate normalization

---

## Design Principles

BOTAS is designed to provide:

- Accurate bacterial RNA-seq alignment
- Native support for circular genomes
- Reproducible gene and operon quantification
- Transparent and interpretable alignment scoring
- Modular architecture for method development and extension
- Integration of alignment and operon-level analyses within a unified workflow

---

## Requirements

### Required

- Python ≥ 3.10
- pysam
- edlib
- biopython

### Optional

- tqdm (progress display)

---

## Documentation

Command-line help is available through:

```bash
botas --help
botas align --help
botas quant --help
botas getOperons --help
```

---

## Citation

If you use BOTAS in your research, please cite:

> Wekesa, C. S. *BOTAS: Bacterial Operon-Aware Transcriptome Alignment System.*

Citation details will be updated following publication.

---

## License

BOTAS is distributed under the MIT License. See the `LICENSE` file for details.

---

## Author

**Clabe Simiyu Wekesa**

GitHub: https://github.com/clabe-wekesa