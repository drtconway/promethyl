import pandas as pd
import pyranges as pr

# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def annotate_methylation(merged_df: pd.DataFrame, ann_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate island-level trio df with gene/promoter information."""
    ann_clean = (
        ann_df
        .drop(columns=["prom_chrom", "prom_start", "prom_end"], errors="ignore")
        .rename(columns={"strand": "gene_strand"})
    )

    sites_pr = pr.PyRanges(
        merged_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    ann_pr = pr.PyRanges(
        ann_clean.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    hits = sites_pr.join(ann_pr, strandedness=False)

    agg = (
        hits.df
        .groupby(["Chromosome", "Start", "End"], as_index=False)
        .agg({
            "gene":       lambda x: ";".join(sorted(set(x.dropna()))),
            "transcript": lambda x: ";".join(sorted(set(x.dropna()))),
            "gene_id":    lambda x: ";".join(sorted(set(x.dropna()))),
            # cpg_island removed — already present from aggregate_to_islands()
        })
    )

    result = (
        merged_df
        .merge(
            agg,
            left_on=["chrom", "start", "end"],
            right_on=["Chromosome", "Start", "End"],
            how="left",
        )
        .drop(columns=["Chromosome", "Start", "End"], errors="ignore")
    )

    for col in ["gene", "transcript", "gene_id"]:
        result[col] = result[col].fillna("")

    return result