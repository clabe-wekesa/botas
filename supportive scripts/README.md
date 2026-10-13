# BOTAS manuscript reproducibility scripts

This directory is a publication-oriented subset of the development scripts used for the BOTAS manuscript. Run from the project root so that paths such as `datasets/`, `BOTAS_references/`, and `results_*` resolve consistently.

The supporting benchmark datasets and reproducibility scripts for BOTAS are archived on Zenodo at https://zenodo.org/records/22101516.

## Exact experimental datasets

The manuscript uses these paired-end SRA runs:

| Sample | SRA |
|---|---|
| B_subtilis_168 | SRR23693988 |
| C_vibrioides_NA1000 | SRR33305313 |
| E_coli_K12 | SRR10059464 |
| M_smegmatis_mc2155 | SRR10377454 |
| P_aeruginosa_PAO1 | SRR10620245 |
| S_aureus_N315 | SRR13251221 |
| V_cholerae_N16961 | SRR37042651 |

Use `00_setup/download_exact_sra.sh` rather than a discovery-based downloader, because exact accession pinning is required for reproducibility. Use `00_setup/prepare_references.py` to retrieve the version-pinned RefSeq assemblies and annotations.

## Table 1: rRNA filtering

1. `01_rrna/build_rrna_benchmark.sh` creates the balanced 100,000-pair truth-aware datasets.
2. Run the three `*_simulated` scripts.
3. Run the three `score_*_simulated.py` scripts.
4. Experimental filtering is reproduced with the three `*_experimental` scripts.

BOTAS rRNA filtering uses `k=18` and `H_min=10` in the current implementation and in the manuscript benchmarks. The reproducibility scripts therefore pass `--rrna-k 18 --rrna-min-hits 10` explicitly so the benchmark does not depend on software-version defaults.

## Table 2: simulated alignment

Use `02_alignment_simulated/run_all_alignment_benchmarks.sh`. It calls `readSimulator_truth_1M.py` and `compare_aligners.py` and produces:

- `master_alignment_metrics.tsv` (junction and non-junction precision/recall/F1)
- `master_resource_usage.tsv` (alignment runtime and peak RAM)
- `master_flagstat.tsv`

The manuscript simulation design is 1,000,000 paired fragments per genome, including 100,000 deliberately generated junction challenge fragments, with a ±10-bp mate-start tolerance.

## Table 3: experimental alignment and junction recovery

1. Run `03_alignment_experimental/run_4_aligners.sh` on the exact seven FASTQ pairs. It produces `results_align/alignment_benchmark_summary.tsv`.
2. Run `03_alignment_experimental/count_junction_fragments.py --bam-root results_align` (or the root containing the final BOTAS BAMs) to count unique XC:Z:JUNCTION-supporting fragments.

Only the unique fragment count is needed for the manuscript table; MAPQ-stratified fields in the development counter are diagnostic and are not part of the final manuscript result.

## Table 4: quantification

Run, in order:

1. `04_quantification/run_botas_quantify_all.sh`
2. `04_quantification/run_featurecounts_all.sh`
3. `04_quantification/run_htseq_all.sh`
4. `python3 04_quantification/summarize_quantification_all.py`

The summary file is `quantification_tool_summary.tsv`.

## Tables 5 and 6: operon inference

1. Build DBTBS truth with `scrape_dbtbs_operons.py` followed by `filter_dbtbs_strict.py`.
2. Build RegulonDB truth with `make_regulondb_ecoli_truth.py` using the downloaded RegulonDB OperonSet table.
3. Prepare Rockhopper annotations with `make_rockhopper_annotations.py` and run Rockhopper. Set `ROCKHOPPER_JAR=/path/to/Rockhopper.jar` before running `run_rockhopper_all.sh`.
4. Run BOTAS with `run_botas_getoperons_all.sh`. This publication copy explicitly passes the manuscript thresholds: same strand, max IGD 70 bp, minimum mean coverage 0.5, minimum coverage ratio 0.25, and minimum score 0.0.
5. `make_table5_curated.py` evaluates predictions against **evaluable annotated protein-coding curated pairs**, matching the manuscript definition and avoiding the older development script's additional coverage-observability filter.
6. Run `collect_botas_operons.py` and then `compare_botas_rockhopper_all.py` for the seven-genome concordance analysis used in Table 6.

## Scripts intentionally excluded

The original development archive also contained temporary cleanup scripts, SGE launchers, old single-genome comparison scripts, duplicate rRNA-build scripts, and MAPQ-fix rerun helpers. Those are not required to reproduce the final manuscript and should not be deposited as the main reproducibility workflow.

## Software versions

Run `00_setup/capture_versions.sh` in the environment used for the final reproduction and deposit the resulting `software_versions.txt` with the outputs.
