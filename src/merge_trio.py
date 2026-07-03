import pandas as pd

# ---------------------------------------------------------------------------
# Merge trio
# ---------------------------------------------------------------------------

def merge_trio(
    proband_df: pd.DataFrame,
    father_df: pd.DataFrame,
    mother_df: pd.DataFrame,
) -> pd.DataFrame:
    """Inner-join the three samples on chrom/start/end.

    Only sites covered at the required depth in all three members are kept.
    """
    coords = ["chrom", "start", "end"]
    merged = (
        proband_df
        .merge(father_df, on=coords, how="inner")
        .merge(mother_df, on=coords, how="inner")
    )
    print(f"  Sites shared across all three samples: {len(merged):,}")
    return merged