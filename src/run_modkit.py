import subprocess
import sys
from pathlib import Path

import pandas as pd

from bedmethyl import (
    read_bedmethyl,
    filter_sites,
    filter_modification,
)

# ---------------------------------------------------------------------------
# modkit
# ---------------------------------------------------------------------------

def run_modkit(
    bam: Path,
    output: Path,
    ref: str=None,
    threads: int=10,
    region: str=None,
    include_bed: str | None=None
) -> None:
    """Run modkit pileup on a single BAM to produce a bedMethyl file."""

    cmd=[
        "modkit", "pileup",
        str(bam),
        str(output),
        "--threads", str(threads),
        "--no-filtering",
        "--motif", "CG", "0",
    ]

    if ref:
        cmd += ["--ref", ref]
    if region:
        cmd += ["--region", region]
    if include_bed:
        cmd += ["--include-bed", include_bed]

    print(f"    $ {' '.join(cmd)}")
    result=subprocess.run(cmd, stderr=subprocess.PIPE)

    if result.returncode != 0:
        raise RuntimeError(f"modkit failed on {bam.name}:\n{result.stderr.decode()}")

# ---------------------------------------------------------------------------
# Per-sample loading
# ---------------------------------------------------------------------------

def ensure_modkit_output(label: str, bam: Path, modkit_dir: Path, ref: str, threads: int, region: str, include_bed) -> Path| None:
    """Run modkit if required and return the bedMethyl output path."""
    modkit_out=modkit_dir/f"{bam.stem}_modkit.bed"
    if modkit_out.exists():
        print(f"  [{label}] modkit output exists → {modkit_out.name}")
    else:
        print(f"  [{label}] running modkit pileup...")
        run_modkit(bam, modkit_out, ref=ref, threads=threads, region=region, include_bed=include_bed)
        print(f"  [{label}] written → {modkit_out.name}")
    if not modkit_out.exists() or modkit_out.stat().st_size == 0:
        print(f"  [{label}] WARNING: no CpG data written for {bam.name} — skipping")
        return None
    return modkit_out


def load_sample_methylation(label: str, bam: Path, modkit_out: Path, min_coverage: int, mod_code: str) -> pd.DataFrame | None:
    """Load, filter and rename a sample's bedMethyl data."""
    df   =read_bedmethyl(str(modkit_out))
    n_raw=len(df)

    if n_raw == 0:
        print(f"  [{label}] WARNING: no CpG sites found in {bam.name} — skipping")
        return None

    df   =filter_sites(df, min_coverage=min_coverage)
    df   =filter_modification(df, mod_code=mod_code)
    print(f"  [{label}] {n_raw:,} → {len(df):,} sites  (cov≥{min_coverage}, mod='{mod_code}')")

    return (
        df[["chrom", "start", "end", "methylation", "coverage", "n_mod", "n_canonical"]]
        .rename(columns={
            "methylation": f"methylation_{label}",
            "coverage":    f"coverage_{label}",
            "n_mod":       f"n_mod_{label}",
            "n_canonical": f"n_canonical_{label}",
        })
    )


def get_sample_methylation(
    label: str,
    bam: Path,
    modkit_dir: Path,
    ref: str,
    threads: int,
    min_coverage: int,
    mod_code: str,
    region: str,
    include_bed=None,
) -> pd.DataFrame | None:
    """Run modkit if needed, load and filter bedMethyl, return per-sample df.

    Columns returned:
        chrom, start, end,
        methylation_<label>, coverage_<label>,
        n_mod_<label>, n_canonical_<label>
    """
    modkit_out=ensure_modkit_output(label=label, bam= bam, modkit_dir= modkit_dir, ref= ref, threads=threads, region= region, include_bed=include_bed)
    if modkit_out is None:
        return None
    df=load_sample_methylation(label =label, bam=bam, modkit_out=modkit_out, min_coverage= min_coverage, mod_code= mod_code)

    return df


