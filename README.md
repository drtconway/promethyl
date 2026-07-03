# promethyl

CpG methylation analysis pipeline for PacBio HiFi long-read sequencing data. Runs `modkit pileup` on haplotagged BAMs (hg38), aggregates methylation to CpG islands, annotates against gene promoters, and tests for differential methylation — either within a trio (proband/father/mother) or across a cohort (outlier detection).

## Requirements

- Conda / mamba
- A reference genome FASTA (hg38), not included in this repo
- `bedtools` on your `PATH` (used when building the promoter annotation reference)

Create the environment:

```bash
conda env create -f environment.yml
conda activate promethyl
```

Key pinned tools/packages: `modkit=0.6.3`, `pandas`, `numpy`, `scipy`, `pyranges`, `pyyaml`.

## One-time setup: build the CpG island / promoter annotation

Before running the pipeline, generate the annotation BED file that both pipeline modes rely on:

```bash
python src/annotate_promoter_cpgs.py
```

This intersects a GENCODE GTF (promoters defined as −2000/+500 bp around each transcript TSS) with a CpG islands BED using `bedtools intersect`, producing `CpGs_with_promoters.bed`.


## Usage

There are two ways to run `src/main.py`: **trio mode** and **cohort mode**.

### Option 1 — Trio mode

Compares a proband against both parents and calls differentially methylated CpG islands.

```bash
python src/main.py \
    --proband proband.bam \
    --father  father.bam \
    --mother  mother.bam \
    --modkit-dir  /data/modkit \
    --annotation  reference/CpGs_with_promoters.bed \
    --ref         /data/genome.fa \
    --output      trio_methylation.tsv \
    --threads 10
```

### Option 2 — Cohort mode

Builds a methylation matrix across many samples and flags outliers by z-score.

```bash
python src/main.py \
    --config sample.yml \
    --modkit-dir /data/modkit/ \
    --annotation reference/CpGs_with_promoters.bed \
    --ref /data/genome.fa \
    --output /output-dir/probands_modkit_cohort_methylation.tsv \
    --include-bed reference/CPGIslandsBED3.bed \
    --threads 10
```

`sample.yml` should list each sample's ID and BAM path:

```yaml
samples:
  - id: sample01
    bam: /data/bams/sample01.bam
  - id: sample02
    bam: /data/bams/sample02.bam
```

### Convenience wrapper

`run_pipeline.sh` wraps cohort mode with repo-relative default paths:

```bash
./run_pipeline.sh sample.yml
```

Edit the `REF` variable near the top of the script to point at your local reference FASTA before running.

`modkit pileup` is skipped automatically for any sample whose bedMethyl output already exists in `--modkit-dir`, so reruns only reprocess new samples.

## Pipeline overview

1. **`run_modkit.py`** — runs `modkit pileup` per sample (skipped if cached output exists)
2. **`bedmethyl.py`** — parses bedMethyl output, filters by coverage and modification code (`m` = 5mC, `h` = 5hmC)
3. **`CpG_meth.py`** — aggregates site-level calls to CpG islands and annotates with overlapping genes/promoters
4. **`merge_trio.py`** — merges proband/father/mother methylation for trio mode
5. **`cohort.py`** — builds the cross-sample matrix and detects per-sample outliers by z-score (cohort mode)
6. **`annotate.py`** — attaches gene/promoter annotation to results
7. **`statistics.py`** — trio-mode significance testing (Fisher's exact / chi-squared) with Benjamini–Hochberg FDR correction

## Key CLI options

| Flag | Description | Default |
|---|---|---|
| `--min-coverage` | Minimum read coverage per CpG site | `10` |
| `--mod-code` | Modification code to test (`m`=5mC, `h`=5hmC) | `m` |
| `--min-delta` | Minimum absolute methylation difference (proband vs. parent mean) | `0.3` |
| `--fdr` | FDR significance threshold (trio mode) | `0.01` |
| `--z-threshold` | Z-score threshold for outlier calling (cohort mode) | `2.0` |
| `--region` | Restrict `modkit pileup` to a genomic region (e.g. `chr1:1000000-1100000`) | — |
| `--include-bed` | Restrict `modkit pileup` to regions in a BED file | — |
| `--threads` | Threads for `modkit pileup` | `10` |

## Output

A TSV of CpG islands with per-sample methylation, coverage, gene/promoter annotation, and (depending on mode) either significance test results or outlier calls.

## Notes

- Outputs, BAMs, BED/TSV files, and the `reference/` directory are git-ignored — regenerate them locally rather than committing.
- `--modkit-dir`, BAMs, and the reference FASTA are expected to live outside the repo.
