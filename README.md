# BOTAS

## Bacterial Operon-aware Transcriptome Alignment Suite

BOTAS is a seed-and-extend RNA-seq alignment and quantification framework designed specifically for bacterial genomes. It combines bacterial read alignment, native circular-genome handling, gene-level quantification, and operon inference in a unified Python command-line application.

Unlike general-purpose RNA-seq aligners designed primarily around eukaryotic splicing, BOTAS targets features of prokaryotic transcriptomes such as dense gene organization, polycistronic transcription, strand specificity, and circular chromosomes and plasmids.

## Key features

- Seed-and-extend alignment optimized for bacterial genomes
- Paired-end and single-end RNA-seq alignment
- Native support for circular chromosomes and plasmids
- BOTAS-native reusable reference indexes
- Optional rRNA filtering using a bundled reference set
- Multiprocessing support for alignment and coverage calculation
- Coordinate-sorted BAM generation and indexing through `pysam`
- Gene-level and operon-level expression quantification
- Operon inference from RNA-seq coverage and genomic organization
- Multi-BAM consensus operon inference
- Docker and Apptainer/Singularity support for reproducible deployment

## Requirements

BOTAS requires Python 3.10 or newer. Its core Python dependencies are:

- `pysam >= 0.22`
- `edlib >= 1.3.9`
- `biopython >= 1.81`

`tqdm >= 4.66` is available as an optional progress-display dependency.

BOTAS does not require external aligners such as Bowtie2, STAR, HISAT2, or BWA for its native alignment workflow.

## Installation

### PyPI

```bash
python -m pip install botas-rnaseq
```

The installed command is:

```bash
botas --help
```

### From source

```bash
git clone https://github.com/clabe-wekesa/botas.git
cd botas
python -m pip install .
```

To include the optional progress display:

```bash
python -m pip install ".[progress]"
```

For development:

```bash
python -m pip install -e ".[dev]"
```

## Command overview

```text
botas
├── index       Build a BOTAS-native reference index
├── align       Align RNA-seq reads to a bacterial reference
├── quantify    Quantify gene or operon expression from BAM files
└── getOperons  Infer bacterial operons from RNA-seq alignments
```

Use `--help` with any command to see all options:

```bash
botas index --help
botas align --help
botas quantify --help
botas getOperons --help
```

## Quick start

### 1. Build a reference index

For a linear reference:

```bash
botas index \
  -r genome.fna \
  -o genome.botas.idx
```

For a circular bacterial reference:

```bash
botas index \
  -r genome.fna \
  --circular \
  --circular-overhang-percent 5 \
  -o genome.circular.botas.idx
```

Individual contigs can instead be marked as circular with `--circular-contigs`.

### 2. Align paired-end reads

```bash
botas -d analysis \
  align \
  -x genome.botas.idx \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -t 8 \
  --pool \
  --sort-bam \
  -o sample.bam
```

With `--sort-bam`, BOTAS also creates a coordinate-sorted BAM and `.bai` index.

For single-end reads, use `--fq`:

```bash
botas -d analysis \
  align \
  -x genome.botas.idx \
  --fq reads.fastq.gz \
  -t 8 \
  --pool \
  --sort-bam \
  -o sample.bam
```

BOTAS can also align directly from a FASTA reference using `-r/--ref`, although a saved BOTAS index is preferable for repeated analyses.

### 3. Quantify gene expression

```bash
botas -d analysis \
  quantify \
  -b analysis/results/sample.sorted.bam \
  -g genome.gff \
  -o sample_counts.tsv
```

The default feature type is `gene`, and the default identifier attribute for genes is `locus_tag`.

BOTAS can also quantify operon features:

```bash
botas -d analysis \
  quantify \
  -b analysis/results/sample.sorted.bam \
  -g operons.gff \
  --feature-type operon \
  --id-attribute ID \
  -o operon_counts.tsv
```

### 4. Infer operons

```bash
botas -d analysis \
  getOperons \
  -b analysis/results/sample.sorted.bam \
  -g genome.gff \
  -t 8 \
  --write-gff
```

Operon inference uses genomic adjacency, strand consistency, and RNA-seq coverage coherence. Multiple BAM files can be analyzed independently or combined using `--consensus`.

## Docker

Official BOTAS container images are published through GitHub Container Registry:

```text
ghcr.io/clabe-wekesa/botas
```

For reproducible analyses, use a specific release tag rather than `latest`.

### Pull the image

```bash
docker pull ghcr.io/clabe-wekesa/botas:0.1.6
```

Test it:

```bash
docker run --rm ghcr.io/clabe-wekesa/botas:0.1.6 --help
```

### Run BOTAS on local files

Mount the current analysis directory at `/work` inside the container:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  ghcr.io/clabe-wekesa/botas:0.1.6 \
  index \
  -r /work/genome.fna \
  -o /work/genome.botas.idx
```

The `--user` option is recommended on Linux so files created by Docker remain owned by the current user.

A paired-end alignment can then be run as:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  ghcr.io/clabe-wekesa/botas:0.1.6 \
  -d /work/analysis \
  align \
  -x /work/genome.botas.idx \
  -1 /work/reads_R1.fastq.gz \
  -2 /work/reads_R2.fastq.gz \
  -t 8 \
  --pool \
  --sort-bam \
  -o sample.bam
```

### Build the image locally

From the BOTAS repository root:

```bash
docker build -t botas:0.1.6 .
docker run --rm botas:0.1.6 --help
```

## Apptainer / Singularity

Apptainer can consume the same public GHCR image, so no separate container definition is required.

### Pull a SIF image

```bash
apptainer pull botas_0.1.6.sif \
  docker://ghcr.io/clabe-wekesa/botas:0.1.6
```

Test it:

```bash
apptainer exec botas_0.1.6.sif botas --help
```

### Run with local data

```bash
apptainer exec \
  --bind "$PWD:/work" \
  botas_0.1.6.sif \
  botas index \
  -r /work/genome.fna \
  -o /work/genome.botas.idx
```

This workflow is particularly useful on HPC systems where Apptainer or Singularity is available but users do not have permission to install system-wide software.

## Container release policy

A GitHub Actions workflow builds and smoke-tests the container whenever a GitHub Release is published. A release tagged, for example, `v0.1.7` publishes:

```text
ghcr.io/clabe-wekesa/botas:0.1.7
ghcr.io/clabe-wekesa/botas:latest
```

`latest` is updated only for non-prerelease releases. For reproducibility, scientific workflows should pin a numbered BOTAS release.

## Architecture

```text
botas/
├── core/       Alignment engine, indexing, pairing, and circular handling
├── io/         FASTQ, BAM, reference, and working-directory handling
├── operons/    Operon inference and consensus logic
├── quantify/   Gene and operon quantification
├── rrna/       rRNA detection and bundled rRNA reference support
├── data/       Packaged runtime data
└── cli/        Command-line interface
```

The alignment engine includes k-mer/minimizer indexing, seed clustering, edit-distance extension using `edlib`, CIGAR construction, pairing logic, MAPQ scoring, and circular-coordinate handling.

## Reproducibility

For reproducible analyses, record the BOTAS version used:

```bash
botas --version
```

For containerized analyses, pin the image tag, for example:

```text
ghcr.io/clabe-wekesa/botas:0.1.6
```

rather than relying on `latest`.

## Citation

If you use BOTAS in your research, please cite:

> Wekesa, C. S. (2026). BOTAS: Bacterial Operon-aware Transcriptome Alignment System. Software manuscript in preparation.

The citation information should be updated when the BOTAS manuscript is published.

## License

BOTAS is distributed under the MIT License. See [LICENSE](LICENSE).

## Author

**Clabe Simiyu Wekesa**  
GitHub: [clabe-wekesa](https://github.com/clabe-wekesa)
