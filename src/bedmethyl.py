from pathlib import Path
import pandas as pd

BEDMETHYL_COLUMNS = [
    "chrom",
    "start",
    "end",
    "mod_code",
    "score",
    "strand",
    "compat_start",
    "compat_end",
    "color",
    "Nvalid_cov",
    "fraction_modified",
    "Nmod",
    "Ncanonical",
    "Nother_mod",
    "Ndelete",
    "Nfail",
    "Ndiff",
    "Nnocall",
]

def read_bedmethyl(filepath):
    """
    Read modkit bedMethyl output into a pandas dataframe.
    """

    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=BEDMETHYL_COLUMNS,
    )

    df = df.rename(
        columns={
            "Nvalid_cov": "coverage",
            "fraction_modified": "methylation",
            "Nmod": "n_mod",
            "Ncanonical": "n_canonical",
        }
    )

    return df[
        [
            "chrom",
            "start",
            "end",
            "mod_code",
            "strand",
            "coverage",
            "methylation",
            "n_mod",
            "n_canonical",
        ]
    ]

def filter_sites(df, min_coverage=10):
    return df[df["coverage"] >= min_coverage].copy()

def filter_modification(df,mod_code="m"):
    return df[df["mod_code"] == mod_code].copy()