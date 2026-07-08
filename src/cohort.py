import numpy as np
import pandas as pd
import pyranges as pr

def bh_adjust(pvals: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    arr = np.asarray(pvals, dtype=float)
    order = np.argsort(arr)
    padj = np.empty(n)

    cummin = 1.0
    for i in range(n - 1, -1, -1):
        cummin = min(cummin, arr[order[i]] * n / (i + 1))
        padj[order[i]] = cummin

    return np.clip(padj, 0.0, 1.0)


def load_config(path: str) -> pd.DataFrame:
    """Load YAML cohort config into a sample metadata dataframe."""
    import yaml
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return pd.DataFrame(cfg["samples"])


def build_cohort_matrix(sample_dfs: dict) -> pd.DataFrame:
    """
    Merge per-sample island methylation into a single wide matrix.
    One row per island, one set of columns per sample.
    """
    coords = ["chrom", "start", "end", "cpg_island"]
    matrix = None

    for sample_id, df in sample_dfs.items():
        # Keep only columns that exist for this sample
        keep = coords + [c for c in df.columns if c.endswith(f"_{sample_id}")]
        sub  = df[[c for c in keep if c in df.columns]]
        if matrix is None:
            matrix = sub
        else:
            matrix = matrix.merge(sub, on=coords, how="inner")

    return matrix


def detect_outliers(
    matrix: pd.DataFrame,
    sample_ids: list[str],
    min_delta: float = 0.20,
    z_threshold: float = 2.0,
) -> pd.DataFrame:
    """
    For each island, compute Z-score of every sample against the
    cohort distribution. Any sample can be flagged as an outlier.
    """
    df = matrix.copy()
    meth_cols = [f"methylation_{s}" for s in sample_ids]

    # Cohort-level stats
    df["cohort_mean"]   = df[meth_cols].mean(axis=1)
    df["cohort_std"]    = df[meth_cols].std(axis=1)
    df["cohort_median"] = df[meth_cols].median(axis=1)

    # Per-sample outlier stats — every sample treated equally
    for sid in sample_ids:
        meth_col = f"methylation_{sid}"
        df[f"delta_{sid}"]   = df[meth_col] - df["cohort_mean"]
        df[f"zscore_{sid}"]  = (
            (df[meth_col] - df["cohort_mean"]) /
            df["cohort_std"].replace(0, np.nan)
        )
        df[f"outlier_{sid}"] = (
            (df[f"zscore_{sid}"].abs() >= z_threshold) &
            (df[f"delta_{sid}"].abs()  >= min_delta)
        )

    # Summary columns
    outlier_cols = [f"outlier_{sid}" for sid in sample_ids]
    df["n_outliers"]      = df[outlier_cols].sum(axis=1)
    df["outlier_samples"] = df[outlier_cols].apply(
        lambda row: ";".join(
            sid for sid, is_out in zip(sample_ids, row) if is_out
        ),
        axis=1,
    )
    df["any_outlier"] = df["n_outliers"] > 0

    return df.sort_values(
        ["n_outliers", "cohort_std"], ascending=[False, False]
    ).reset_index(drop=True)


def to_long_format(df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """
    Reshape a wide cohort dataframe (one row per island, per-sample columns
    suffixed with "_<sample_id>") into long format: one row per
    (island, sample), with per-sample metrics as plain-named columns.

    Locus/annotation/cohort-level columns (chrom, start, end, cpg_island,
    gene, cohort_mean, ...) are not sample-specific and are repeated across
    each island's sample rows.
    """
    per_sample_prefixes = [
        col[: -(len(sid) + 1)]
        for sid in sample_ids
        for col in df.columns
        if col.endswith(f"_{sid}")
    ]
    per_sample_prefixes = sorted(set(per_sample_prefixes))

    shared_cols = [
        col for col in df.columns
        if not any(col == f"{prefix}_{sid}" for prefix in per_sample_prefixes for sid in sample_ids)
    ]

    frames = []
    for sid in sample_ids:
        sample_cols = {f"{prefix}_{sid}": prefix for prefix in per_sample_prefixes if f"{prefix}_{sid}" in df.columns}
        sub = df[shared_cols + list(sample_cols.keys())].rename(columns=sample_cols)
        sub.insert(len(shared_cols), "sample", sid)
        frames.append(sub)

    long_df = pd.concat(frames, ignore_index=True)

    if "outlier_samples" in long_df.columns:
        long_df = long_df.drop(columns=["outlier_samples"])

    sort_cols = [c for c in ["chrom", "start", "end", "sample"] if c in long_df.columns]
    return long_df.sort_values(sort_cols).reset_index(drop=True)


def compute_cohort_statistics(
    matrix: pd.DataFrame,
    sample_ids: list[str],
    min_delta: float = 0.20,
    z_threshold: float = 2.0,
) -> pd.DataFrame:
    """
    For every sample at every island, run a one-vs-rest chi2 test.
    Every sample is treated as a potential outlier.
    """
    from scipy.stats import chi2_contingency

    df = matrix.copy()

    def chi2_p(row, n_mod_a, n_can_a, n_mod_b, n_can_b):
        table = [
            [row[n_mod_a], row[n_can_a]],
            [row[n_mod_b], row[n_can_b]],
        ]
        if any(sum(r) == 0 for r in table) or any(sum(c) == 0 for c in zip(*table)):
            return 1.0
        _, p, _, _ = chi2_contingency(table, correction=False)
        return p

    for sid in sample_ids:
        rest = [s for s in sample_ids if s != sid]

        # Pool all other samples
        df[f"_n_mod_rest"]  = df[[f"n_mod_{s}"       for s in rest]].sum(axis=1)
        df[f"_n_can_rest"]  = df[[f"n_canonical_{s}" for s in rest]].sum(axis=1)

        df[f"pval_{sid}"] = df.apply(
            chi2_p, axis=1,
            n_mod_a=f"n_mod_{sid}",
            n_can_a=f"n_canonical_{sid}",
            n_mod_b="_n_mod_rest",
            n_can_b="_n_can_rest",
        )
        df[f"padj_{sid}"] = bh_adjust(df[f"pval_{sid}"])
        df = df.drop(columns=["_n_mod_rest", "_n_can_rest"])

    return df