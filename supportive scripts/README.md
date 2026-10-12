# BOTAS manuscript reproducibility scripts

This package contains the scripts required to reproduce the analyses reported in the BOTAS manuscript. The executable scripts do **not** contain fixed sample names, absolute paths, or manuscript-specific directory names. Input datasets are supplied through command-line arguments or tab-separated manifest files, so the same scripts can be used with the manuscript datasets or with other bacterial RNA-Seq datasets.

The manuscript-specific SRA accessions and RefSeq assemblies are stored separately under `config/`. These are data manifests, not hard-coded logic in the analysis scripts.

## 1. Requirements

The complete workflow uses the following programs, depending on the analysis being reproduced:

- Python 3.9 or newer
- BOTAS
- samtools
- SRA Toolkit (`fasterq-dump`)
- curl
- ART/ART Modern (`art_modern`)
- bedtools
- Bowtie2
- BWA-MEM2
- minimap2
- SortMeRNA
- RiboDetector (`ribodetector_cpu`)
- Subread/featureCounts
- HTSeq
- Java and Rockhopper

Python scripts that inspect BAM files require `pysam`.

Record the installed versions before starting:

```bash
bash 00_setup/capture_versions.sh software_versions.txt
```

## 2. Manifest files

Most scripts operate on TSV manifests. Paths may be absolute or relative to the directory from which the command is executed.

### 2.1 Experimental sample manifest

The principal manifest has five columns:

```text
sample    r1    r2    reference    gff
```

An example is provided as `config/samples.template.tsv`.

`sample` is an arbitrary unique sample label. `r1` and `r2` are paired-end FASTQ files. `reference` is the reference FASTA used for alignment. `gff` is the corresponding GFF3 annotation.

### 2.2 BAM manifest

Downstream quantification and operon scripts use:

```text
sample    bam    gff    reference
```

An example is provided as `config/bams.template.tsv`.

### 2.3 Manuscript-specific accession manifests

The exact SRA runs are in:

```text
config/manuscript_sra.tsv
```

The exact RefSeq assemblies are in:

```text
config/manuscript_references.tsv
```

These files can be replaced with manifests for any other samples without editing the scripts.

## 3. Download the manuscript data

Create a working directory first:

```bash
mkdir -p work
```

### 3.1 Download the seven paired-end SRA datasets

```bash
python3 00_setup/download_sra.py \
    --manifest config/manuscript_sra.tsv \
    --outdir work/fastq \
    --threads 16
```

The output files are named from the `sample` column of the manifest, for example:

```text
work/fastq/B_subtilis_168_1.fastq.gz
work/fastq/B_subtilis_168_2.fastq.gz
```

The same downloader can be used for any paired-end SRA runs by supplying another TSV with columns `sample` and `sra`.

### 3.2 Download the version-pinned reference genomes and GFF3 files

```bash
python3 00_setup/prepare_references.py \
    --manifest config/manuscript_references.tsv \
    --out work/references
```

Each sample receives its own directory containing the prepared FASTA and GFF3 files. Checksums and reference metadata are written to `work/references/references.tsv` and `work/references/replicons.tsv`.

### 3.3 Build the analysis sample manifest

```bash
python3 00_setup/build_samples_manifest.py \
    --samples config/manuscript_sra.tsv \
    --fastq-dir work/fastq \
    --reference-dir work/references \
    --output work/samples.tsv
```

Inspect the file before running benchmarks:

```bash
column -t -s $'\t' work/samples.tsv
```

## 4. Table 1: rRNA filtering

The manuscript uses BOTAS with `k = 18` and `H_min = 10`. Both parameters are passed explicitly below so the reproduction does not depend on a future software default.

### 4.1 Build the truth-aware simulated rRNA benchmark

```bash
bash 01_rrna/build_rrna_benchmark.sh \
    --manifest work/samples.tsv \
    --outdir work/rrna_sim \
    --rrna-pairs 50000 \
    --non-rrna-pairs 50000 \
    --read-length 150 \
    --insert-mean 300 \
    --insert-sd 50 \
    --threads 16
```

This produces, for every sample, 50,000 canonical-rRNA read pairs and 50,000 non-rRNA read pairs together with a truth table.

Build a manifest for those simulated FASTQ files:

```bash
python3 01_rrna/build_rrna_sim_manifest.py \
    --source-manifest work/samples.tsv \
    --benchmark-dir work/rrna_sim \
    --output work/rrna_sim_samples.tsv
```

### 4.2 Run BOTAS on the simulated benchmark

```bash
python3 01_rrna/run_rrna_tools.py \
    --tool botas \
    --manifest work/rrna_sim_samples.tsv \
    --outdir work/rrna_tools \
    --threads 16 \
    --rrna-k 18 \
    --rrna-min-hits 10
```

### 4.3 Run SortMeRNA

Supply the FASTA used as the SortMeRNA rRNA reference database:

```bash
python3 01_rrna/run_rrna_tools.py \
    --tool sortmerna \
    --manifest work/rrna_sim_samples.tsv \
    --outdir work/rrna_tools \
    --rrna-ref /path/to/sortmerna_rRNA_database.fasta \
    --threads 16
```

### 4.4 Run RiboDetector

```bash
python3 01_rrna/run_rrna_tools.py \
    --tool ribodetector \
    --manifest work/rrna_sim_samples.tsv \
    --outdir work/rrna_tools \
    --read-length 150 \
    --threads 16
```

Each invocation writes an output manifest such as:

```text
work/rrna_tools/botas_outputs.tsv
work/rrna_tools/sortmerna_outputs.tsv
work/rrna_tools/ribodetector_outputs.tsv
```

### 4.5 Score the simulated rRNA benchmark

BOTAS:

```bash
python3 01_rrna/score_rrna.py \
    --truth-manifest work/rrna_sim_samples.tsv \
    --predictions-manifest work/rrna_tools/botas_outputs.tsv \
    --output work/rrna_tools/botas_metrics.tsv
```

SortMeRNA:

```bash
python3 01_rrna/score_rrna.py \
    --truth-manifest work/rrna_sim_samples.tsv \
    --predictions-manifest work/rrna_tools/sortmerna_outputs.tsv \
    --output work/rrna_tools/sortmerna_metrics.tsv
```

RiboDetector:

```bash
python3 01_rrna/score_rrna.py \
    --truth-manifest work/rrna_sim_samples.tsv \
    --predictions-manifest work/rrna_tools/ribodetector_outputs.tsv \
    --output work/rrna_tools/ribodetector_metrics.tsv
```

### 4.6 Experimental rRNA filtering

Run the same three tools using `work/samples.tsv` instead of the simulated manifest. For BOTAS:

```bash
python3 01_rrna/run_rrna_tools.py \
    --tool botas \
    --manifest work/samples.tsv \
    --outdir work/rrna_experimental \
    --threads 16 \
    --rrna-k 18 \
    --rrna-min-hits 10
```

The SortMeRNA and RiboDetector commands are identical to those above except for the input manifest and output directory. Experimental data do not have a truth table, so these runs are summarized by the number or percentage of read pairs classified as rRNA-like rather than by precision, recall, or F1. For each tool, summarize the run using its output manifest. For example, for BOTAS:

```bash
python3 01_rrna/summarize_rrna_experimental.py \
    --samples work/samples.tsv \
    --predictions-manifest work/rrna_experimental/botas_outputs.tsv \
    --output work/rrna_experimental/botas_summary.tsv
```

Repeat with `sortmerna_outputs.tsv` and `ribodetector_outputs.tsv` for the other tools.

## 5. Table 2: truth-aware simulated alignment benchmark

The simulation/alignment driver accepts any directory containing FASTA files. It recursively searches the supplied directory; sample names are derived from FASTA filenames.

For the manuscript design:

```bash
bash 02_alignment_simulated/run_all_alignment_benchmarks.sh \
    --reference-dir work/references \
    --outdir work/alignment_simulated \
    --threads 16 \
    --total-pairs 1000000 \
    --junction-pairs 100000 \
    --iterations 10 \
    --read-length 150 \
    --seed 12345 \
    --tolerance 10 \
    --overhang-percent 5
```

The principal outputs are:

```text
work/alignment_simulated/master_alignment_metrics.tsv
work/alignment_simulated/master_resource_usage.tsv
work/alignment_simulated/master_flagstat.tsv
```

### 5.1 `readSimulator_truth_1M.py`: generate truth-aware paired-end data

This script is the foundation of the simulated alignment benchmark. It accepts any FASTA containing one or more circular replicons, generates paired-end reads with classic ART, and writes both FASTQ files and the mate-level truth table used for accuracy evaluation.

For the manuscript-style design:

```bash
python3 02_alignment_simulated/readSimulator_truth_1M.py \
    --input work/references/E_coli_K12/reference.fna \
    --outdir work/simulated/E_coli_K12 \
    --iterations 10 \
    --readlen 150 \
    --total-pairs 1000000 \
    --junction-pairs 100000 \
    --opts "-m 300 -s 30" \
    --seed 12345
```

If `art_illumina` is not available as `art` on `PATH`, provide the executable explicitly with `--simulator_path /path/to/art_illumina`.

Important arguments:

| Argument | Meaning |
|---|---|
| `--input` | FASTA containing the circular replicon(s) to simulate |
| `--simulator_path` | ART executable; default `art` |
| `--outdir` | Output directory |
| `--iterations` | Number of random whole-replicon rotations used for background simulation |
| `--readlen` | Read length in bp |
| `--total-pairs` | Exact number of read pairs in the final dataset |
| `--junction-pairs` | Number of truth-confirmed origin-spanning challenge pairs retained in the final dataset |
| `--opts` | Additional paired-end ART options; the manuscript uses `"-m 300 -s 30"` |
| `--seed` | Random seed for reproducibility |

For an input FASTA named `reference.fna`, principal outputs include:

```text
reference_1.fastq
reference_2.fastq
reference_truth.tsv
```

The truth TSV contains the expected reference sequence, 1-based start coordinate, strand and origin-crossing status for each mate, together with the reference length. This file is the authoritative truth input for the evaluation scripts below.

```bash
python3 02_alignment_simulated/readSimulator_truth_1M.py -h
```

### 5.2 `compare_aligners.py`: compare several BAMs against the same truth

Use this script to reproduce the manuscript comparison of BOTAS with other aligners. The same truth table is applied to every BAM, and metrics are reported separately for `ALL`, `JUNCTION`, and `NON_JUNCTION` reads.

```bash
python3 02_alignment_simulated/compare_aligners.py \
    --truth work/simulated/E_coli_K12/reference_truth.tsv \
    --bam BOTAS work/bams/E_coli_K12.botas.sorted.bam \
    --bam Bowtie2 work/bams/E_coli_K12.bowtie2.sorted.bam \
    --bam BWA-MEM2 work/bams/E_coli_K12.bwa_mem2.sorted.bam \
    --bam minimap2 work/bams/E_coli_K12.minimap2.sorted.bam \
    --tolerance 10 \
    --tsv work/results/E_coli_K12_alignment_metrics.tsv
```

Each `--bam` takes two values: a label and a BAM path. Repeat the option for as many aligners as required.

A mapped mate is counted as correct only when its reference sequence and strand match truth and its reported 1-based start coordinate lies within `--tolerance` bp of the true start. The manuscript uses 10 bp.

The categories are defined from simulation truth:

- `ALL`: every simulated mate;
- `JUNCTION`: truth-labelled origin-spanning mates;
- `NON_JUNCTION`: all remaining truth-labelled mates.

`XC:Z:JUNCTION` is not used to decide whether a simulated read belongs to the junction category. If present, the number of such tags is reported separately as a BOTAS diagnostic.

The optional `--tsv` output contains `aligner`, `category`, `total`, `correct`, `wrong`, `unmapped`, `precision`, `recall`, `f1`, and `junction_tags`.

```bash
python3 02_alignment_simulated/compare_aligners.py -h
```

### 5.3 `evaluate_pe_bam.py`: evaluate one BAM independently

Use `evaluate_pe_bam.py` for single-BAM validation and debugging. Unlike `compare_aligners.py`, it accepts truth supplied as SAM or BAM in addition to the BOTAS simulator TSV.

```bash
python3 02_alignment_simulated/evaluate_pe_bam.py \
    --bam work/bams/E_coli_K12.botas.sorted.bam \
    --truth work/simulated/E_coli_K12/reference_truth.tsv \
    --out work/results/E_coli_K12_botas_evaluation.tsv \
    --tolerance 10
```

Supported truth inputs are:

```text
.tsv   BOTAS truth table produced by readSimulator_truth_1M.py
.sam   SAM containing the known paired-end alignments
.bam   BAM containing the known paired-end alignments
```

Important arguments:

| Argument | Meaning |
|---|---|
| `--bam` | BAM being evaluated |
| `--truth` | Truth TSV, SAM or BAM |
| `--out` | Output metrics TSV |
| `--tolerance` | Maximum permitted start-coordinate deviation in bp; default 10 |

The output contains counts of correct, wrong and unmapped reads together with precision, recall and F1. A summary is also printed to the terminal.

Use `compare_aligners.py` when several aligners should be compared in one command with separate junction/non-junction summaries. Use `evaluate_pe_bam.py` for single-BAM validation.

```bash
python3 02_alignment_simulated/evaluate_pe_bam.py -h
```

## 6. Table 3: experimental alignment benchmark

Run all four aligners on any samples defined in `samples.tsv`:

```bash
python3 03_alignment_experimental/run_4_aligners.py \
    --manifest work/samples.tsv \
    --outdir work/alignment_experimental \
    --threads 24 \
    --circular all \
    --botas-overhang-percent 5
```

The script writes:

```text
work/alignment_experimental/alignment_benchmark_summary.tsv
```

This file gives the BAM path for every sample/aligner combination.

### 6.1 Create a BOTAS BAM manifest

```bash
python3 00_setup/select_bams.py \
    --samples work/samples.tsv \
    --alignment-summary work/alignment_experimental/alignment_benchmark_summary.tsv \
    --aligner BOTAS \
    --output work/botas_bams.tsv
```

### 6.2 Count experimental circular-junction fragments

```bash
python3 03_alignment_experimental/count_junction_fragments.py \
    --manifest work/botas_bams.tsv \
    --output work/junction_fragment_summary.tsv \
    --details work/junction_fragment_details.tsv \
    --min-mapq 20
```

The manuscript uses the unique fragment count carrying `XC:Z:JUNCTION`; the MAPQ-stratified field is diagnostic.

## 7. Table 4: gene quantification

Use the BOTAS BAM manifest produced above:

```bash
python3 04_quantification/run_quantification.py \
    --manifest work/botas_bams.tsv \
    --outdir work/quantification \
    --threads 24 \
    --feature-type gene \
    --id-attribute locus_tag \
    --min-overlap 1
```

This runs BOTAS quantify, featureCounts and HTSeq using the same BAM/GFF pairs.

Summarize assignment counts:

```bash
python3 04_quantification/summarize_quantification_all.py \
    --manifest work/botas_bams.tsv \
    --results-dir work/quantification \
    --output work/quantification_tool_summary.tsv
```

## 8. Tables 5 and 6: operon inference

### 8.1 Run BOTAS getOperons

The manuscript thresholds are supplied explicitly:

```bash
python3 05_operons/run_botas_getoperons.py \
    --manifest work/botas_bams.tsv \
    --outdir work/botas_operons \
    --threads 24 \
    --gene-id gene \
    --max-igd 70 \
    --min-coverage 0.5 \
    --min-cov-ratio 0.25 \
    --min-score 0.0
```

The script writes `work/botas_operons/operon_outputs.tsv`, which maps each sample to its BOTAS operon TSV.

### 8.2 DBTBS curated truth set

First scrape/export the DBTBS material using `scrape_dbtbs_operons.py` according to the source format used in the manuscript. The scraper itself exposes its options with:

```bash
python3 05_operons/scrape_dbtbs_operons.py -h
```

Then apply the strict experimental-evidence filter:

```bash
python3 05_operons/filter_dbtbs_strict.py \
    --input dbtbs_operons_experimental.tsv \
    --operons-out dbtbs_operons_strict.tsv \
    --pairs-out dbtbs_operon_strict_truth_pairs.tsv
```

### 8.3 RegulonDB curated truth set

After obtaining the RegulonDB `OperonSet` table, run:

```bash
python3 05_operons/make_regulondb_ecoli_truth.py \
    --operonset /path/to/OperonSet.tsv \
    --gff /path/to/E_coli_annotation.gff3 \
    --operons-out regulondb_operons_clean.tsv \
    --pairs-out regulondb_truth_pairs.tsv \
    --unmapped-out regulondb_unmapped_members.tsv
```

No fixed E. coli path is embedded in the script.

### 8.4 Prepare Rockhopper annotations

```bash
python3 05_operons/make_rockhopper_annotations.py \
    --manifest work/samples.tsv \
    --outdir work/rockhopper_references
```

Rockhopper requires its own reference/annotation representation. Build the Rockhopper run manifest automatically:

```bash
python3 05_operons/build_rockhopper_manifest.py \
    --samples work/samples.tsv \
    --reference-dir work/rockhopper_references \
    --output work/rockhopper_samples.tsv
```

The resulting manifest has columns `sample`, `r1`, `r2`, and `rockhopper_reference`. Then run:

```bash
python3 05_operons/run_rockhopper.py \
    --manifest work/rockhopper_samples.tsv \
    --outdir work/rockhopper \
    --jar /path/to/Rockhopper.jar \
    --threads 24 \
    --java-mem 16g
```

### 8.5 Curated comparison for Table 5

Create a TSV with one row per curated organism and these columns:

```text
sample    reference_name    gff    truth_pairs    botas_operons    rockhopper_operons
```

Then run:

```bash
python3 05_operons/make_table5_curated.py \
    --manifest work/table5_inputs.tsv \
    --output work/table5_curated_operon.tsv
```

The comparison is restricted to curated adjacent pairs whose genes are represented as protein-coding genes in the supplied annotation. Predictions absent from the curated database are reported as `predictions_not_represented`; they are not automatically treated as biological false positives.

### 8.6 BOTAS–Rockhopper concordance for Table 6

Create a TSV containing:

```text
sample    gff    botas_operons    rockhopper_operons
```

Then run:

```bash
python3 05_operons/compare_botas_rockhopper_all.py \
    --manifest work/table6_inputs.tsv \
    --output work/botas_rockhopper_concordance.tsv
```

The script compares adjacent protein-coding gene pairs while preserving original adjacency; it does not reconnect genes across intervening non-protein-coding features.

## 9. Running one sample instead of the complete manuscript panel

The scripts do not require seven organisms. Make a manifest containing a single row and run exactly the same commands. For example:

```text
sample    r1    r2    reference    gff
my_sample /data/R1.fastq.gz /data/R2.fastq.gz /data/genome.fna /data/genome.gff3
```

This is also the recommended way to test the workflow before launching all manuscript datasets.

## 10. Re-running and output safety

Most batch runners skip outputs that already exist. Use `--force` where available when a benchmark parameter, BOTAS version, reference, or input FASTQ has changed. Do not combine outputs from runs produced with different parameters.

For publication reproduction, keep the following with the final output archive:

- `software_versions.txt`
- the manifests actually used
- `work/references/references.tsv` and checksum fields
- benchmark result TSVs
- the exact BOTAS version/commit or package version

## 11. Checking script interfaces

Every Python executable provides command-line help:

```bash
python3 path/to/script.py -h
```

The two shell drivers provide:

```bash
bash 01_rrna/build_rrna_benchmark.sh --help
bash 02_alignment_simulated/run_all_alignment_benchmarks.sh --help
```

The key design rule of this package is that **file locations and sample identities live in manifests or command-line arguments, not inside the analysis code**.
