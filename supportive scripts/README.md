# BOTAS manuscript supplementary data and benchmarking scripts

This repository contains the simulated datasets, evaluation scripts, and benchmarking workflows used in the manuscript describing **BOTAS**, a bacterial RNA-seq analysis framework.

The archive is intended to complement the BOTAS software repository by providing the data and analysis procedures required to reproduce the benchmarking analyses reported in the manuscript.

## BOTAS software

BOTAS itself is distributed separately from this supplementary archive.

The analyses in this archive require an installed version of BOTAS together with the external programs used for individual benchmarking comparisons.

---

## Archive structure

```text
.
├── botas_manuscript_data/
│   └── simulated/
│       ├── circular_wrap/
│       │   ├── ecoli/
│       │   ├── rphaseoli/
│       │   └── rphaseoli_se/
│       ├── standard_pe/
│       │   ├── ecoli/
│       │   └── rphaseoli/
│       └── standard_se/
│           └── rphaseoli/
│
├── botas_manuscript_scripts/
│   ├── benchmark_operon_links_unordered.py
│   ├── compare_botas_featurecounts.py
│   ├── evaluate_pe_bam.py
│   ├── evaluate_se_bam.py
│   ├── extract_regulondb_operon_links.py
│   ├── plot_operon_benchmark.py
│   ├── prepare_rphaseoli_circular_truth.py
│   ├── run_final_ecoli_circular_pe.sh
│   ├── run_final_ecoli_standard_pe_updated.sh
│   ├── run_final_rphaseoli_circular_pe.sh
│   ├── run_final_rphaseoli_standard_pe.sh
│   ├── run_final_rphaseoli_standard_se.sh
│   ├── run_quantification_benchmark.sh
│   ├── simulate_rphaseoli_circular_wrap.py
│   └── simulate_rphaseoli_circular_wrap_se.py
│
├── references/
├── index/
└── results/
```

The `references/`, `index/`, and `results/` directories are used by the benchmarking scripts. Reference genomes and external datasets should be obtained from the sources described in the manuscript and placed at the locations indicated below.

---

# 1. Requirements

The scripts were developed for a Linux environment.

The main software requirements are:

- BOTAS
- Python 3
- pysam
- pandas
- NumPy
- SciPy
- scikit-learn
- matplotlib
- pymongo/BSON support
- samtools
- Bowtie2
- BWA-MEM2
- Minimap2
- HISAT2
- STAR
- featureCounts

GNU `/usr/bin/time` is used for runtime and peak-memory measurements.

The benchmarking shell scripts use six threads by default. The number of threads can be changed using the `THREADS` environment variable, for example:

```bash
THREADS=8 ./botas_manuscript_scripts/run_final_ecoli_standard_pe_updated.sh
```

---

# 2. Reference genomes

The benchmarking workflows expect the following directory structure:

```text
references/
├── ecoli/
│   ├── genome.fa
│   └── annotation.gff
└── rphaseoli/
    └── genome.fa
```

The *Escherichia coli* and *Rhizobium phaseoli* reference sequences used for the manuscript analyses are described in the Methods section of the manuscript.

The *R. phaseoli* reference corresponds to strain R650, which contains one chromosome and four plasmids.

---

# 3. Simulated datasets

The simulated reads used for alignment benchmarking are included under:

```text
botas_manuscript_data/simulated/
```

Three principal benchmark types are represented:

```text
standard_pe/
circular_wrap/
standard_se/
```

`standard_pe` contains conventional paired-end simulated reads distributed across the reference genome.

`circular_wrap` contains reads or fragments specifically designed to cross the end/origin boundary of circular bacterial replicons.

`standard_se` contains the single-end benchmark dataset.

The corresponding truth tables contain the expected reference sequence and genomic alignment positions used to evaluate mapping accuracy.

---

# 4. Paired-end alignment evaluation

Paired-end BAM files are evaluated with:

```text
botas_manuscript_scripts/evaluate_pe_bam.py
```

The evaluator compares primary alignments with the expected reference sequence and genomic start position in the truth table.

An alignment is considered correct when:

1. it maps to the expected reference sequence; and
2. its start position is within the specified positional tolerance.

The default tolerance is 10 bp.

Secondary and supplementary alignments are excluded.

The script reports:

- total read pairs;
- total reads;
- mapped reads;
- correctly mapped reads;
- incorrectly mapped reads;
- unmapped reads;
- precision;
- recall; and
- F1 score.

Example:

```bash
python botas_manuscript_scripts/evaluate_pe_bam.py \
    --bam results/alignment/ecoli_standard_pe/bowtie2.bam \
    --truth botas_manuscript_data/simulated/standard_pe/ecoli/ecoli_20x_truth.tsv \
    --out results/metrics/ecoli_standard_pe/bowtie2.metrics.tsv \
    --tolerance 10
```

---

# 5. Single-end alignment evaluation

Single-end alignments are evaluated using:

```text
botas_manuscript_scripts/evaluate_se_bam.py
```

Example:

```bash
python botas_manuscript_scripts/evaluate_se_bam.py \
    --bam <alignment.bam> \
    --truth botas_manuscript_data/simulated/standard_se/rphaseoli/rphaseoli_se_truth.tsv \
    --out <metrics.tsv> \
    --tolerance 10
```

The supplied single-end truth tables use 0-based genomic coordinates, corresponding directly to `pysam.AlignmentSegment.reference_start`.

---

# 6. Standard paired-end alignment benchmarks

## E. coli

Run:

```bash
./botas_manuscript_scripts/run_final_ecoli_standard_pe_updated.sh
```

The workflow compares:

- BOTAS
- Bowtie2
- BWA-MEM2
- Minimap2
- HISAT2
- STAR

using the same simulated paired-end dataset.

The summary is written to:

```text
results/tables/ecoli_standard_pe_summary.tsv
```

## R. phaseoli

Run:

```bash
./botas_manuscript_scripts/run_final_rphaseoli_standard_pe.sh
```

The corresponding summary is written to:

```text
results/tables/rphaseoli_standard_pe_summary.tsv
```

For each aligner, the scripts record runtime, peak RAM, precision, recall, and F1 score.

---

# 7. Circular-boundary paired-end benchmarks

Circular-wrap datasets test the ability of an aligner to recover reads originating from fragments spanning the end/origin boundary of circular bacterial replicons.

## E. coli

Run:

```bash
./botas_manuscript_scripts/run_final_ecoli_circular_pe.sh
```

## R. phaseoli

Run:

```bash
./botas_manuscript_scripts/run_final_rphaseoli_circular_pe.sh
```

BOTAS is run with circular-genome handling enabled:

```bash
botas align --circular ...
```

The other aligners are evaluated against the same reads and truth positions.

---

# 8. R. phaseoli circular truth-table conversion

The raw *R. phaseoli* circular-wrap truth table stores simulated positions as 0-based coordinates:

```text
rphaseoli_wrap_truth.tsv
```

The paired-end evaluator uses the 1-based `start1_1` and `start1_2` convention used by the standard paired-end truth tables.

The evaluation truth table is therefore generated with:

```bash
python botas_manuscript_scripts/prepare_rphaseoli_circular_truth.py \
    --input botas_manuscript_data/simulated/circular_wrap/rphaseoli/rphaseoli_wrap_truth.tsv \
    --output botas_manuscript_data/simulated/circular_wrap/rphaseoli/rphaseoli_wrap_truth_eval.tsv
```

The resulting file is:

```text
botas_manuscript_data/simulated/circular_wrap/rphaseoli/rphaseoli_wrap_truth_eval.tsv
```

This conversion adds one nucleotide to the raw 0-based positions without otherwise changing the simulated truth information.

---

# 9. R. phaseoli single-end benchmark

Run:

```bash
./botas_manuscript_scripts/run_final_rphaseoli_standard_se.sh
```

The workflow compares BOTAS, Bowtie2, BWA-MEM2, Minimap2, HISAT2, and STAR.

The BOTAS single-end benchmark uses:

```text
--step 15
```

as used in the manuscript analysis.

---

# 10. Circular-wrap simulation scripts

Scripts used for generating *R. phaseoli* circular-boundary reads are provided under:

```text
botas_manuscript_scripts/
```

The single-end circular-wrap dataset can be generated using:

```bash
python botas_manuscript_scripts/simulate_rphaseoli_circular_wrap_se.py
```

The default simulation parameters are:

```text
reads per replicon: 10,000
read length:         150 bp
random seed:         42
```

Every generated read is forced to cross the end/origin boundary of its source replicon.

A paired-end circular simulation script is also supplied as:

```text
simulate_rphaseoli_circular_wrap.py
```

The archived simulated FASTQ and truth files used for the reported benchmark are supplied directly under `botas_manuscript_data/`.

---

# 11. Gene quantification benchmark

BOTAS gene-level quantification was compared with featureCounts using the same paired-end experimental *E. coli* alignment file and genome annotation.

Run:

```bash
BAM=/path/to/ecoli.sorted.bam \
GFF=references/ecoli/annotation.gff \
THREADS=6 \
./botas_manuscript_scripts/run_quantification_benchmark.sh
```

BOTAS is run using:

```bash
botas quantify \
    -b <BAM> \
    -g <GFF> \
    --gff-gene-attribute ID
```

featureCounts is run in paired-fragment counting mode using:

```text
-p --countReadPairs
```

Runtime and peak memory are measured for both programs.

Outputs are written under:

```text
results/quantification/
```

---

# 12. BOTAS versus featureCounts agreement

After running the quantification benchmark, calculate gene-level agreement using:

```bash
python botas_manuscript_scripts/compare_botas_featurecounts.py
```

The script calculates:

- number of genes compared;
- number of genes with identical counts;
- Pearson correlation;
- Spearman correlation;
- R²;
- mean absolute error (MAE);
- root mean squared error (RMSE); and
- maximum absolute count difference.

The default outputs are:

```text
results/quantification/botas_featurecounts_comparison.tsv
results/quantification/botas_featurecounts_gene_counts.tsv
```

The second file contains the gene-by-gene BOTAS and featureCounts values underlying the reported statistics.

---

# 13. Operon benchmark against RegulonDB

BOTAS operon predictions were evaluated against experimentally curated transcription-unit information from RegulonDB.

The benchmark consists of three stages.

## 13.1 Extract RegulonDB transcription units

Run:

```bash
python botas_manuscript_scripts/extract_regulondb_operon_links.py \
    --bson /path/to/operonDatamart.bson
```

This produces:

```text
results/operons/regulondb/
├── regulondb_adjacent_gene_links.tsv
└── regulondb_transcription_units.tsv
```

Adjacent genes within each RegulonDB transcription unit are treated as reference operon links.

The RegulonDB release used for the manuscript analysis should be matched to the version reported in the manuscript.

## 13.2 Compare BOTAS with RegulonDB

A gene-name-to-gene-ID mapping is required because RegulonDB gene names and the annotation identifiers used by BOTAS are not necessarily identical.

The mapping file must contain:

```text
gene_name	gene_id
```

Run:

```bash
python botas_manuscript_scripts/benchmark_operon_links_unordered.py \
    --gene-map references/ecoli/regulondb_gene_map.tsv \
    --regulondb-links results/operons/regulondb/regulondb_adjacent_gene_links.tsv \
    --botas-operons results/operons/botas_operons.tsv
```

Gene-pair direction is ignored during comparison. Thus:

```text
geneA-geneB
```

and:

```text
geneB-geneA
```

are treated as the same adjacent-gene link. This prevents reverse-strand transcription units from being penalized solely because their genes are represented in the opposite order.

The benchmark reports:

- mapped reference links;
- predicted links;
- true positives;
- false positives;
- false negatives;
- unmapped reference links;
- precision;
- recall; and
- F1 score.

Detailed TP, FP, and FN link tables are also retained.

## 13.3 Plot the operon benchmark

Run:

```bash
python botas_manuscript_scripts/plot_operon_benchmark.py
```

Figures are written to:

```text
results/operons/figures/
```

in both PNG and PDF formats.

---

# 14. Experimental RNA-seq data

Raw experimental RNA-seq reads are not redistributed in this archive because they are publicly available from the Sequence Read Archive.

The manuscript analyses include experimental datasets for *E. coli* and *R. phaseoli*. Their SRA accession numbers and reference genome accessions are provided in the manuscript.

Users should retrieve the original reads directly from the corresponding public repository before reproducing the experimental-data analyses.

---

# 15. Output directories

The supplied workflows create results under:

```text
results/
├── alignment/
├── metrics/
├── tables/
├── quantification/
├── experimental/
└── operons/
```

Large intermediate BAM files and aligner indices do not need to be archived because they can be regenerated from the supplied data and reference sequences.

---

# 16. Reproducibility notes

The simulated FASTQ files and truth tables used for the manuscript benchmarks are included directly in this archive.

Random seeds are fixed where applicable.

Alignment accuracy is evaluated against known simulated genomic positions rather than inferred from alignment statistics alone.

Runtime and peak-memory measurements use GNU `/usr/bin/time`.

Benchmark scripts accept the `THREADS` environment variable where appropriate.

Reference genomes, external databases, and publicly archived experimental RNA-seq datasets should be obtained from their original repositories using the accessions and versions reported in the manuscript.

---

# 17. Citation

If you use these supplementary benchmarking materials, please cite the associated BOTAS manuscript and the archived Zenodo dataset.

Wekesa, C., Kelvin, K., Muoma, J., & Mithöfer, A. (2026). *BOTAS: Supplementary Benchmarking Data and Scripts for Bacterial RNA-seq Analysis* (Version 1.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.22101516

The BOTAS manuscript will also be cited once its publication details are available.

---

# 18. License

The supplementary benchmarking data and scripts in this archive are released under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

Copyright © 2026 Clabe Wekesa.

Under the CC BY 4.0 license, these materials may be shared and adapted provided appropriate credit is given to the creator(s), a link to the license is provided, and any changes are indicated.

License: https://creativecommons.org/licenses/by/4.0/

Zenodo DOI: https://doi.org/10.5281/zenodo.22101516

---