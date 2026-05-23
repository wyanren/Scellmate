#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


# =========================================================
# args
# =========================================================
parser = argparse.ArgumentParser(
    description="Per-eMGE per-species one-sided Fisher enrichment filtering for candidate linkages."
)
parser.add_argument(
    "-i", "--input",
    default="candidate_linkages_0.7score.tsv",
    help="Input candidate linkage matrix (default: candidate_linkages_0.7score.tsv)"
)
parser.add_argument(
    "-s", "--species",
    default="final.reformat.add_cluster.edit_cont.annotate.species.txt",
    help="SAG species annotation file (default: final.reformat.add_cluster.edit_cont.annotate.species.txt)"
)
parser.add_argument(
    "--alpha",
    type=float,
    default=0.05,
    help="P-value cutoff for one-sided Fisher test (default: 0.05)"
)
parser.add_argument(
    "--min_count",
    type=int,
    default=1,
    help="Minimum number of positive calls required within a species for testing (default: 1)"
)

# default ON
parser.add_argument(
    "--remove_singleton",
    dest="remove_singleton",
    action="store_true",
    default=True,
    help="Remove species groups with a_this_pos == 1 even if Fisher p < alpha (default: ON)"
)
parser.add_argument(
    "--keep_singleton",
    dest="remove_singleton",
    action="store_false",
    help="Keep singleton groups if Fisher p < alpha"
)

args = parser.parse_args()

input_file = Path(args.input)
species_file = Path(args.species)
alpha = args.alpha
min_count = args.min_count
remove_singleton = args.remove_singleton

if not input_file.exists():
    raise FileNotFoundError(f"Input file not found: {input_file}")
if not species_file.exists():
    raise FileNotFoundError(f"Species file not found: {species_file}")


# =========================================================
# load species annotation
# =========================================================
sp = pd.read_csv(
    species_file,
    sep=r"\s+",
    header=None,
    names=["cluster", "sag", "genus", "species"]
)

sp = sp.drop_duplicates("sag").set_index("sag")["species"].astype(str)

species_to_short = {
    "Bacillus_subtilis": "B_subtilis",
    "Escherichia_coli": "E_coli",
    "Klebsiella_pneumoniae": "K_pneumoniae",
    "Staphylococcus_aureus": "S_aureus",
}


# =========================================================
# helper: infer true species from eMGE column name
# =========================================================
def infer_emge_species(colname: str):
    for full_name, short_name in species_to_short.items():
        if short_name in colname:
            return full_name
    return None


# =========================================================
# load candidate linkage matrix
# =========================================================
df = pd.read_csv(input_file, sep="\t", index_col=0)

# keep only SAGs with species annotation
common_sags = df.index.intersection(sp.index)
df = df.loc[common_sags].copy()
sp = sp.loc[common_sags].copy()

# make sure all values are 0/1 integers
df = df.fillna(0).astype(int)

filtered = df.copy()
test_records = []
removal_records = []

species_list = sorted(sp.unique())


# =========================================================
# per-column per-species fisher test
# =========================================================
for col in df.columns:
    y = df[col].to_numpy()
    expected_species = infer_emge_species(col)

    for species_name in species_list:
        mask_this = (sp.values == species_name)
        mask_other = ~mask_this

        a = int(y[mask_this].sum())   # this species, linkage=1
        b = int(mask_this.sum() - a)  # this species, linkage=0
        c = int(y[mask_other].sum())  # other species, linkage=1
        d = int(mask_other.sum() - c) # other species, linkage=0

        tested = True
        odds_ratio = np.nan
        pvalue = np.nan
        fisher_passed = False
        removed = False
        removal_reason = None
        final_keep = False
        reason = "tested"

        # skip trivial cases
        if a < min_count:
            tested = False
            reason = f"skip_a_lt_min_count({min_count})"
        elif (a + c) == 0:
            tested = False
            reason = "skip_no_positive_in_column"
        elif mask_this.sum() == 0 or mask_other.sum() == 0:
            tested = False
            reason = "skip_empty_group"
        else:
            odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            fisher_passed = bool(pvalue < alpha)

            if not fisher_passed:
                removed = (a > 0)
                removal_reason = "fisher_failed" if removed else None
                final_keep = False
                reason = "tested_fisher_failed"
            else:
                if remove_singleton and a == 1:
                    removed = True
                    removal_reason = "singleton_removed"
                    final_keep = False
                    reason = "tested_fisher_passed_but_singleton_removed"
                else:
                    removed = False
                    removal_reason = None
                    final_keep = True
                    reason = "tested_passed"

        # if should remove, zero out all linkages for this species in this column
        if tested and removed and a > 0:
            row_idx = sp.index[mask_this]

            # classify removed linkages by ground-truth consistency
            if expected_species is not None and species_name == expected_species:
                removed_consistent = a
                removed_inconsistent = 0
                removed_group_type = "species_consistent"
            else:
                removed_consistent = 0
                removed_inconsistent = a
                removed_group_type = "species_inconsistent"

            removal_records.append({
                "eMGE": col,
                "expected_species_from_eMGE": expected_species,
                "species": species_name,
                "a_this_pos": a,
                "removed_linkages": a,
                "removed_group_type": removed_group_type,
                "removed_species_consistent_linkages": removed_consistent,
                "removed_species_inconsistent_linkages": removed_inconsistent,
                "removal_reason": removal_reason,
                "odds_ratio": odds_ratio,
                "pvalue_one_sided_greater": pvalue,
                "alpha": alpha,
                "remove_singleton": remove_singleton
            })

            filtered.loc[row_idx, col] = 0

        test_records.append({
            "eMGE": col,
            "expected_species_from_eMGE": expected_species,
            "species": species_name,
            "species_size": int(mask_this.sum()),
            "other_size": int(mask_other.sum()),
            "a_this_pos": a,
            "b_this_neg": b,
            "c_other_pos": c,
            "d_other_neg": d,
            "odds_ratio": odds_ratio,
            "pvalue_one_sided_greater": pvalue,
            "alpha": alpha,
            "fisher_passed": fisher_passed,
            "remove_singleton": remove_singleton,
            "singleton_removed": bool(tested and fisher_passed and remove_singleton and a == 1),
            "final_keep": final_keep,
            "removed": removed,
            "removal_reason": removal_reason,
            "tested": tested,
            "reason": reason
        })


# =========================================================
# species-consistent check after filtering
# =========================================================
consistency_records = []

total_retained = int(filtered.to_numpy().sum())
consistent_total = 0
inconsistent_total = 0
unknown_total = 0

for col in filtered.columns:
    expected_species = infer_emge_species(col)
    col_pos = filtered[col] == 1

    if expected_species is None:
        for sag_id in filtered.index[col_pos]:
            consistency_records.append({
                "SAG": sag_id,
                "SAG_species": sp.loc[sag_id],
                "eMGE": col,
                "expected_species_from_eMGE": None,
                "species_consistent": np.nan
            })
            unknown_total += 1
        continue

    for sag_id in filtered.index[col_pos]:
        sag_species = sp.loc[sag_id]
        is_consistent = (sag_species == expected_species)

        consistency_records.append({
            "SAG": sag_id,
            "SAG_species": sag_species,
            "eMGE": col,
            "expected_species_from_eMGE": expected_species,
            "species_consistent": int(is_consistent)
        })

        if is_consistent:
            consistent_total += 1
        else:
            inconsistent_total += 1

consistency_df = pd.DataFrame(consistency_records)

denom = consistent_total + inconsistent_total
summary_df = pd.DataFrame([{
    "total_retained_linkages": total_retained,
    "species_consistent_linkages": consistent_total,
    "species_inconsistent_linkages": inconsistent_total,
    "species_unknown_linkages": unknown_total,
    "species_consistent_fraction":
        (consistent_total / denom) if denom > 0 else np.nan
}])


# =========================================================
# summary of removed eMGE-species groups
# =========================================================
removal_df = pd.DataFrame(removal_records)

if removal_df.shape[0] > 0:
    removed_group_count = int(removal_df.shape[0])
    removed_total_linkages = int(removal_df["removed_linkages"].sum())
    removed_consistent_total = int(removal_df["removed_species_consistent_linkages"].sum())
    removed_inconsistent_total = int(removal_df["removed_species_inconsistent_linkages"].sum())

    removed_group_consistent = int((removal_df["removed_group_type"] == "species_consistent").sum())
    removed_group_inconsistent = int((removal_df["removed_group_type"] == "species_inconsistent").sum())

    removed_by_fisher_groups = int((removal_df["removal_reason"] == "fisher_failed").sum())
    removed_by_singleton_groups = int((removal_df["removal_reason"] == "singleton_removed").sum())

    removed_by_fisher_linkages = int(removal_df.loc[
        removal_df["removal_reason"] == "fisher_failed", "removed_linkages"
    ].sum())
    removed_by_singleton_linkages = int(removal_df.loc[
        removal_df["removal_reason"] == "singleton_removed", "removed_linkages"
    ].sum())
else:
    removed_group_count = 0
    removed_total_linkages = 0
    removed_consistent_total = 0
    removed_inconsistent_total = 0
    removed_group_consistent = 0
    removed_group_inconsistent = 0
    removed_by_fisher_groups = 0
    removed_by_singleton_groups = 0
    removed_by_fisher_linkages = 0
    removed_by_singleton_linkages = 0


# =========================================================
# cleaned refinement result table for Supplementary Data
# =========================================================
test_df = pd.DataFrame(test_records)


def summarize_refinement_result(row):
    a = int(row["a_this_pos"])
    fisher_passed = bool(row["fisher_passed"])
    singleton_removed = bool(row["singleton_removed"])
    final_keep = bool(row["final_keep"])
    tested = bool(row["tested"])

    if final_keep:
        return "Pass"

    if a == 0:
        return "No self candidate linkage"

    if (not tested) and a > 0:
        return "Fail (not tested)"

    if fisher_passed and singleton_removed:
        return "Fail (Fisher passed; singleton)"

    if (not fisher_passed) and a == 1:
        return "Fail (Fisher failed; singleton)"

    if not fisher_passed:
        return "Fail (Fisher failed)"

    return "Fail"


refinement_result_df = test_df.copy()

refinement_result_df["Result"] = refinement_result_df.apply(
    summarize_refinement_result,
    axis=1
)

refinement_result_df = refinement_result_df.rename(columns={
    "species": "Focal species",
    "species_size": "SAG number in focal species",
    "a_this_pos": "Self positive linkages",
    "b_this_neg": "Self negative linkages",
    "c_other_pos": "Other positive linkages",
    "d_other_neg": "Other negative linkages",
    "pvalue_one_sided_greater": "P value",
    "fisher_passed": "Fisher passed",
    "singleton_removed": "Singleton removed",
    "final_keep": "Final kept"
})

# Keep only eMGE-species groups with at least one self candidate linkage.
# These are the candidate linkage groups that were actually refined.
refinement_result_df = refinement_result_df[
    refinement_result_df["Self positive linkages"] > 0
].copy()

refinement_result_df = refinement_result_df[[
    "eMGE",
    "Focal species",
    "SAG number in focal species",
    "Self positive linkages",
    "Self negative linkages",
    "Other positive linkages",
    "Other negative linkages",
    "P value",
    "Fisher passed",
    "Singleton removed",
    "Final kept",
    "Result"
]]


# =========================================================
# write outputs
# =========================================================
prefix = input_file.name.replace(".tsv", "")

filtered_file = input_file.with_name(f"{prefix}.fisher_filtered.tsv")
test_file = input_file.with_name(f"{prefix}.fisher_test_results.tsv")
refinement_result_file = input_file.with_name(f"{prefix}.refinement_result.tsv")
removal_file = input_file.with_name(f"{prefix}.fisher_removed_emge_species.tsv")
consistency_file = input_file.with_name(f"{prefix}.fisher_species_consistency.tsv")
summary_file = input_file.with_name(f"{prefix}.fisher_species_consistency_summary.tsv")

filtered.to_csv(filtered_file, sep="\t")
test_df.to_csv(test_file, sep="\t", index=False)
refinement_result_df.to_csv(refinement_result_file, sep="\t", index=False)
removal_df.to_csv(removal_file, sep="\t", index=False)
consistency_df.to_csv(consistency_file, sep="\t", index=False)
summary_df.to_csv(summary_file, sep="\t", index=False)


# =========================================================
# final echo
# =========================================================
print(f"[out] filtered matrix: {filtered_file}")
print(f"[out] fisher test results: {test_file}")
print(f"[out] cleaned refinement result: {refinement_result_file}")
print(f"[out] removed eMGE-species groups: {removal_file}")
print(f"[out] species consistency detail: {consistency_file}")
print(f"[out] species consistency summary: {summary_file}")

print("\n===== PARAMETERS =====")
print(f"alpha = {alpha}")
print(f"min_count = {min_count}")
print(f"remove_singleton = {remove_singleton}")

print("\n===== RETAINED SUMMARY =====")
print(summary_df.to_string(index=False))

print("\n===== REMOVED eMGE-SPECIES SUMMARY =====")
print(f"removed_eMGE_species_groups = {removed_group_count}")
print(f"removed_total_linkages      = {removed_total_linkages}")
print(f"removed_species_consistent_linkages   = {removed_consistent_total}")
print(f"removed_species_inconsistent_linkages = {removed_inconsistent_total}")
print(f"removed_species_consistent_groups     = {removed_group_consistent}")
print(f"removed_species_inconsistent_groups   = {removed_group_inconsistent}")
print(f"removed_by_fisher_groups              = {removed_by_fisher_groups}")
print(f"removed_by_singleton_groups           = {removed_by_singleton_groups}")
print(f"removed_by_fisher_linkages            = {removed_by_fisher_linkages}")
print(f"removed_by_singleton_linkages         = {removed_by_singleton_linkages}")

