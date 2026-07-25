# BOTAS

## Bacterial Operon-aware Transcriptome Alignment System

BOTAS is a seed-and-extend RNA-seq alignment system designed specifically for bacterial genomes. It integrates circular genome support, operon-aware alignment logic, strand-aware quantification, and operon inference within a unified Python framework.

Unlike general-purpose RNA-seq aligners optimized for eukaryotic splicing, BOTAS is tailored for prokaryotic transcriptomics, where operon structure, dense gene organization, and circular chromosomes require specialized handling.

---

## Key Features

- Seed-and-extend alignment engine optimized for bacterial genomes  
- Native circular chromosome and plasmid support  
- Paired-end and single-end read alignment  
- Operon-aware alignment logic  
- Strand-aware gene quantification  
- Gene-level and operon-level expression quantification  
- Operon inference from aligned reads  
- Optional rRNA filtering  
- Parallel execution for scalable performance  
- Modular and extensible Python architecture  

---

## Why BOTAS?

Most RNA-seq aligners were developed for eukaryotic transcriptomes and assume splicing, large introns, and linear chromosomes. Bacterial transcriptomes differ fundamentally:

- Genes are densely packed  
- Transcripts frequently span multiple genes (operons)  
- Genomes are often circular  
- Strand specificity is biologically critical  

BOTAS addresses these characteristics directly by integrating circular-aware coordinate handling, operon-consistent pairing logic, insert-size validation in circular space, and integrated gene and operon quantification.

---

## Installation

### Install from PyPI

pip install botas

### Install from source

git clone https://github.com/cwekesa/botas  
cd botas  
pip install .

### Development installation

pip install -e .[dev]

---

## Basic Usage

### Paired-end Alignment

botas align \
  -r reference.fasta \
  -1 reads_R1.fastq \
  -2 reads_R2.fastq \
  -o aligned.bam

### Single-end Alignment

botas align \
  -r reference.fasta \
  -U reads.fastq \
  -o aligned.bam

---

## Gene Quantification

botas quant \
  -b aligned.bam \
  -g genes.gff \
  -o gene_counts.tsv

Supported features include:

- Strand-specific counting (unstranded, forward, reverse)  
- MAPQ filtering  
- Multi-mapper handling (ignore, unique-only, fractional counting)  
- TPM and RPKM calculation  
- Multi-BAM count matrix generation  

---

## Operon Inference

botas getOperons \
  -b aligned.bam \
  -g genes.gff \
  -o operons.tsv

Operon inference integrates:

- Intergenic distance  
- Strand consistency  
- Coverage similarity  
- Paired-end support  
- Consensus-based merging  

---

## Architecture

botas/
├── core/          # Alignment engine and k-mer index
├── io/            # FASTQ, BAM, and reference handling
├── operons/       # Operon inference logic
├── quantify/      # Gene and operon quantification
├── rrna/          # rRNA detection and filtering
└── cli/           # Command-line interface

The alignment engine implements:

- K-mer indexing  
- Seed clustering  
- Edit-distance extension (edlib-based)  
- CIGAR reconstruction  
- MAPQ scoring  
- Circular coordinate normalization  

---

## Design Principles

BOTAS is built with the following goals:

- High-precision bacterial RNA-seq alignment  
- Explicit support for circular genomes  
- Transparent and interpretable scoring  
- Modular architecture for research extensibility  
- Reproducible gene and operon quantification  
- Integration of alignment and operon-level biology  

---

## Requirements

- Python ≥ 3.10  
- pysam  
- edlib  
- biopython  

Optional:

- tqdm (progress display)  
- pytest (development and testing)  

---

## Citation

If you use BOTAS in your research, please cite:

Wekesa, C. S. (2026). BOTAS: Bacterial Operon-aware Transcriptome Alignment System. (Software manuscript in preparation.)

---

## License

MIT License

---

## Author

Clabe Simiyu Wekesa  
GitHub: https://github.com/cwekesa
