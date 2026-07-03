import pandas as pd
import subprocess

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
# gunzip gencode.v49.annotation.gtf.gz

gtf = "/Users/Tim/Documents/Scripts/triomethyl/Reference/gencode.v49.annotation.gtf"

cols = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
]

df = pd.read_csv(gtf, sep="\t", comment="#", names=cols)

# keep transcripts only
df = df[df["feature"] == "transcript"]

def parse_attr(attr, key):
    for x in attr.split(";"):
        x = x.strip()
        if x.startswith(key):
            return x.split('"')[1]
    return None

df["gene_id"] = df["attribute"].apply(lambda x: parse_attr(x, "gene_id"))
df["gene_name"] = df["attribute"].apply(lambda x: parse_attr(x, "gene_name"))
df["transcript_id"] = df["attribute"].apply(lambda x: parse_attr(x, "transcript_id"))

# compute TSS (strand-aware)
df["tss"] = df.apply(
    lambda r: r["start"] if r["strand"] == "+" else r["end"],
    axis=1
)

# promoter window: -2000 to +500
def promoter(row):
    if row["strand"] == "+":
        start = y if (y:= row["tss"] - 2000) > 0 else 0
        end = row["tss"] + 500
    else:
        start = y if (y:= row["tss"] - 500) > 0 else 0
        end = row["tss"] + 2000
    return pd.Series([start, end])

df[["p_start", "p_end"]] = df.apply(promoter, axis=1)

# clean BED output
out = df[[
    "chrom",
    "p_start",
    "p_end",
    "gene_name",
    "transcript_id",
    "strand",
    "gene_id"
]]

# BED is 0-based
out.to_csv(
    "/Users/Tim/Documents/Scripts/triomethyl/Reference/GENCODEV49_promoters_2kb_plus500.bed",
    sep="\t",
    header=False,
    index=False
)

promoters_bed = "/Users/Tim/Documents/Scripts/triomethyl/Reference/GENCODEV49_promoters_2kb_plus500.bed"
cpg_bed = "/Users/Tim/Documents/Scripts/triomethyl/Reference/CPGIslands.bed"  # your input CpG file
output = "/Users/Tim/Documents/Scripts/triomethyl/Reference/CpGs_in_promoters.bed"

cmd = [
    "bedtools", "intersect",
    "-a", cpg_bed,
    "-b", promoters_bed,
    "-wa", "-wb",
]

with open("/Users/Tim/Documents/Scripts/triomethyl/Reference/CpGs_with_promoters.bed", "w") as out:
    subprocess.run(cmd, stdout=out, check=True)