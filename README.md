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

### Parallel cohort mode (Nextflow)

For cohort mode, `main.nf` parallelizes phase 1 (`modkit pileup` + island aggregation) across all samples, then runs phase 2 (cross-sample outlier analysis) once all samples finish:

```bash
nextflow run main.nf --samples sample.yml --outdir results
```

Every pipeline setting other than `--outdir` and `--threads` (per-sample `modkit pileup` threads, still a CLI/`params` option) is read from top-level keys in `sample.yml` instead of `--flag`s:

```yaml
samples:
  - id: sample01
    bam: /data/bams/sample01.bam
  - id: sample02
    bam: /data/bams/sample02.bam

annotation:   reference/CpGs_with_promoters.bed   # required
ref:          /data/genome.fa                     # optional
include_bed:  reference/CPGIslandsBED3.bed         # optional
region:       chr1:1000000-1100000                 # optional
min_coverage: 10
mod_code:     m
min_delta:    0.3
z_threshold:  2.0
```

`annotation` is required in the YAML; the rest fall back to the same defaults as `main.py` (shown above) when omitted.

Each sample runs as its own `MODKIT` process (so samples run concurrently, limited only by the executor's available resources), writing `results/modkit/<sample>_modkit.bed`. Once every sample's process completes, a single `COHORT` process merges their island-level outputs and writes `results/cohort_methylation.tsv`.

`nextflow.config` defines a `standard` profile (local executor, the default) and a `slurm` profile (`-profile slurm`) for running on a Slurm cluster.

Nextflow always uses a pre-built conda environment for process execution — it does not solve `promethyl_environment.yml` on the fly. By default it looks for the env at `envs/promethyl` under `$PROMETHYL_HOME` (falling back to the repo root if `PROMETHYL_HOME` isn't set), so build it there once:

```bash
conda env create -f promethyl_environment.yml -p "${PROMETHYL_HOME:-.}/envs/promethyl"
```

Override the location with `--conda_env /some/other/env` if you'd rather keep it elsewhere.

### Shared install / running batches from elsewhere

To install promethyl once in a shared location and run batches from other directories (e.g. per-cohort working directories on an HPC), rather than checking out the repo per batch:

1. **Check out the repo once** at a shared, read-only-by-convention path, e.g. `/shared/tools/promethyl`.

2. **Set `PROMETHYL_HOME` and pre-build the conda environment there** (recommended on clusters where compute nodes lack channel/internet access):

   ```bash
   export PROMETHYL_HOME=/shared/tools/promethyl
   conda env create -f "$PROMETHYL_HOME/promethyl_environment.yml" -p "$PROMETHYL_HOME/envs/promethyl"
   ```

   `nextflow.config` reads `$PROMETHYL_HOME` at run time (via `System.getenv`), so as long as it's exported before `nextflow run`, the env at `$PROMETHYL_HOME/envs/promethyl` is picked up automatically — no `--conda_env` flag needed.

3. **Use the `bin/promethyl` wrapper** so batch directories don't need to know the full path to `main.nf`:

   ```bash
   export PATH="$PROMETHYL_HOME/bin:$PATH"

   cd /data/cohorts/my-batch
   promethyl --samples cohort.yaml -profile slurm
   ```

   `promethyl` just runs `nextflow run "$PROMETHYL_HOME/main.nf" "$@"`; Nextflow's `work/` directory and outputs still land in the current directory (the batch dir), only the pipeline code and env are shared. Put the `export PROMETHYL_HOME=...` and `export PATH=...` lines in a module file or shared shell profile so users don't set them by hand each time.

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

Trio mode writes one row per CpG island (wide form), with per-sample methylation/coverage columns and significance test results.

Cohort mode (both `main.py --config` and `main.nf`) writes long-form output: one row per (CpG island, sample), with a `sample` column and per-sample metrics (`methylation`, `coverage`, `n_mod`, `n_canonical`, `delta`, `zscore`, `outlier`, `pval`, `padj`) as plain columns. Locus and cohort-level columns (`chrom`, `start`, `end`, `cpg_island`, gene/promoter annotation, `cohort_mean`, `cohort_std`, `cohort_median`, `n_outliers`, `any_outlier`) repeat across each island's sample rows.

## Notes

- Outputs, BAMs, BED/TSV files, and the `reference/` directory are git-ignored — regenerate them locally rather than committing.
- `--modkit-dir`, BAMs, and the reference FASTA are expected to live outside the repo.
