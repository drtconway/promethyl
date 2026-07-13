#!/usr/bin/env python3
"""
run_cohort.py - Phase 2 of the cohort pipeline: outlier analysis across samples.

Collects the per-sample island-methylation TSVs written by run_sample.py
(phase 1), builds the cross-sample matrix, detects per-sample outliers by
z-score, annotates with gene/promoter info, and writes the final cohort TSV.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from CpG_meth import load_promoter_annotations
from annotate import annotate_methylation
from cohort import build_cohort_matrix, detect_outliers, compute_cohort_statistics, to_long_format


def parse_args():
    p = argparse.ArgumentParser(
        description="Cohort outlier analysis across per-sample island methylation TSVs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--samples", required=True, nargs="+",
                   help="id=path pairs for each sample's phase-1 TSV, e.g. sample01=sample01.tsv")
    p.add_argument("--sample-meta", nargs="+", default=None,
                   help="id=tissue=sex triples, e.g. sample01=WB=female. When given, cohort "
                        "statistics for a sample are computed only against other samples of the "
                        "same tissue (and, on chrX, the same sex too). Omit to compare every "
                        "sample against every other sample regardless of tissue/sex.")
    p.add_argument("--annotation", required=True, help="CpGs_with_promoters.bed")
    p.add_argument("--output", required=True, help="Output cohort TSV")

    p.add_argument("--min-delta", type=float, default=0.3, help="Minimum absolute methylation difference")
    p.add_argument("--z-threshold", type=float, default=2.0, help="Z-score threshold for outlier calling")

    return p.parse_args()


def main():
    args = parse_args()

    annotation = Path(args.annotation)
    output = Path(args.output)

    if not annotation.exists():
        sys.exit(
            f"ERROR: annotation file not found: {annotation}\n"
            f"       Run annotate_promoter_cpgs.py first to generate it."
        )

    output.parent.mkdir(parents=True, exist_ok=True)

    ann_df = load_promoter_annotations(str(annotation))
    print(f"  {len(ann_df):,} CpG island-promoter records loaded")

    print(f"=== Cohort mode: {len(args.samples)} samples ===")
    sample_dfs = {}
    for pair in args.samples:
        if "=" not in pair:
            sys.exit(f"ERROR: --samples entries must be id=path, got: {pair}")
        sid, path = pair.split("=", 1)
        path = Path(path)
        if not path.exists():
            sys.exit(f"ERROR: [{sid}] phase-1 TSV not found: {path}")
        sample_dfs[sid] = pd.read_csv(path, sep="\t")
        print(f"  [{sid}] loaded {len(sample_dfs[sid]):,} islands from {path.name}")

    if len(sample_dfs) < 2:
        sys.exit("ERROR: fewer than 2 samples with data — cannot run cohort analysis")

    all_ids = list(sample_dfs.keys())

    sample_meta = None
    if args.sample_meta:
        sample_meta = {}
        for triple in args.sample_meta:
            parts = triple.split("=")
            if len(parts) != 3:
                sys.exit(f"ERROR: --sample-meta entries must be id=tissue=sex, got: {triple}")
            sid, tissue, sex = parts
            sample_meta[sid] = {"tissue": tissue or None, "sex": sex or None}
        missing = [sid for sid in all_ids if sid not in sample_meta]
        if missing:
            print(f"  WARNING: no tissue/sex metadata for {len(missing)} sample(s): {', '.join(missing)}"
                  f" — they will be excluded from every comparison group.")

    matrix = build_cohort_matrix(sample_dfs)
    print(f"\n  Islands covered across all samples: {len(matrix):,}")

    result = compute_cohort_statistics(matrix, all_ids, sample_meta=sample_meta,
                                       min_delta=args.min_delta,
                                       z_threshold=args.z_threshold)
    result = detect_outliers(result, all_ids, sample_meta=sample_meta,
                             min_delta=args.min_delta,
                             z_threshold=args.z_threshold)
    result = annotate_methylation(result, ann_df)
    result = to_long_format(result, all_ids)

    result.to_csv(output, sep="\t", index=False)
    print(f"\n✓ Done!")
    print(f"  Sites tested : {len(result):,}")
    print(f"  Output       : {output}")


if __name__ == "__main__":
    main()
