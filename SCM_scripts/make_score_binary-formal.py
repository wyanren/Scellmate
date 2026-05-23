#!/usr/bin/env python3

import re
import math
import argparse
from pathlib import Path
import pandas as pd


# =========================
# defaults
# =========================
DEFAULT_INPUT_FILE = "all_mock_samples.breadth_absolute.matrix.tsv"
DEFAULT_OUTPUT_FILE = "all_mock_samples.score0.3.binary.tsv"
A0 = 1000.0
DEFAULT_SCORE_CUTOFF = 0.7


# =========================
# helper
# =========================
len_re = re.compile(r"_length_(\d+)(?:_|$)")

def extract_length(colname: str):
    m = len_re.search(str(colname))
    return int(m.group(1)) if m else None


def calc_score(abs_bp: float, length: int, a0: float = 1000.0) -> float:
    if pd.isna(abs_bp) or abs_bp <= 0 or length <= 0:
        return 0.0
    rel = abs_bp / length
    return rel * math.log(1.0 + abs_bp / a0)   # ln instead of log10


# =========================
# arguments
# =========================
parser = argparse.ArgumentParser(
    description="Convert breadth absolute matrix to binary linkage matrix using score threshold."
)
parser.add_argument(
    "-i", "--input",
    default=DEFAULT_INPUT_FILE,
    help=f"Input TSV file (default: {DEFAULT_INPUT_FILE})"
)
parser.add_argument(
    "-o", "--output",
    default=DEFAULT_OUTPUT_FILE,
    help=f"Output TSV file (default: {DEFAULT_OUTPUT_FILE})"
)
parser.add_argument(
    "-s", "--score",
    type=float,
    default=DEFAULT_SCORE_CUTOFF,
    help=f"Score cutoff (default: {DEFAULT_SCORE_CUTOFF})"
)

args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_FILE = args.output
SCORE_CUTOFF = args.score


# =========================
# load matrix
# =========================
infile = Path(INPUT_FILE)
if not infile.exists():
    raise FileNotFoundError(f"Input file not found: {infile}")

df = pd.read_csv(infile, sep="\t", index_col=0)

print(f"[load] matrix shape: {df.shape[0]} rows x {df.shape[1]} MGEs")

# parse lengths from column names
length_map = {}
bad_cols = []
for col in df.columns:
    L = extract_length(col)
    if L is None:
        bad_cols.append(col)
    else:
        length_map[col] = L

if bad_cols:
    print("[warn] Could not parse length for these columns:")
    for c in bad_cols:
        print(f"  {c}")
    print("[warn] These columns will be filled with 0 in output.")

# =========================
# calculate binary linkage
# =========================
out = pd.DataFrame(index=df.index, columns=df.columns, dtype="int8")

for col in df.columns:
    L = length_map.get(col, None)
    if L is None:
        out[col] = 0
        continue

    scores = df[col].apply(lambda x: calc_score(x, L, A0))
    out[col] = (scores >= SCORE_CUTOFF).astype("int8")

# =========================
# summary
# =========================
total_linkages = int(out.to_numpy().sum())
rows_with_linkage = int((out.sum(axis=1) > 0).sum())
mges_with_linkage = int((out.sum(axis=0) > 0).sum())

print(f"[summary] score cutoff = {SCORE_CUTOFF}")
print(f"[summary] total linkage cells = {total_linkages}")
print(f"[summary] rows with >=1 linkage = {rows_with_linkage}")
print(f"[summary] MGEs with >=1 linkage = {mges_with_linkage}")

# =========================
# write output
# =========================
out.to_csv(OUTPUT_FILE, sep="\t")
print(f"[out] wrote: {OUTPUT_FILE}")

# =========================
# final echo
# =========================
print("\n===== FINAL REPORT =====")
print(f"TOTAL_LINKAGES={total_linkages}")
print(f"OUTPUT_FILE={OUTPUT_FILE}")

