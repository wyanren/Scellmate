#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import random
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pick 4 genus-based 3-species genome sets, with at least one Archaea genus."
    )
    parser.add_argument(
        "--taxonomy_tsv",
        default="/mnt/md0/DB/GTDB/database_GTDB_related_full_v220/taxonomy_r220.processed_with_headers.tsv",
        help="GTDB taxonomy table with headers"
    )
    parser.add_argument(
        "--rep_dir",
        default="/mnt/md0/DB/GTDB/database_GTDB_related_full_v220/gtdb_rep/",
        help="Directory containing GTDB representative genome fasta files (*.fna.gz)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed"
    )
    parser.add_argument(
        "--output_tsv",
        default="selected_4_genus_triplets.no_placeholder.tsv",
        help="Output TSV file"
    )
    return parser.parse_args()


def genome_file_exists(rep_dir: Path, genome_id: str) -> bool:
    return (rep_dir / f"{genome_id}_genomic.fna.gz").exists()


def is_placeholder_species(species: str) -> bool:
    """
    Return True for placeholder species names like:
      BM004_sp012974315
      Thermofilum_sp020832965
      WAOL01_sp027023255

    Keep proper named species like:
      Thermofilum_adornatum
      Spiroplasma_A_gladiatoris
    """
    if pd.isna(species):
        return True

    species = str(species).strip()

    # Common GTDB placeholder pattern: *_sp + digits
    if re.search(r"_sp\d+$", species):
        return True

    return False


def load_taxonomy(taxonomy_tsv: Path, rep_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(taxonomy_tsv, sep="\t", dtype=str)

    required_cols = {"Genome_ID", "Species", "Genus", "Domain"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in taxonomy file: {sorted(missing)}")

    df = df.dropna(subset=["Genome_ID", "Species", "Genus", "Domain"]).copy()
    df = df[df["Species"].str.strip() != ""].copy()
    df = df[df["Genus"].str.strip() != ""].copy()

    # Remove placeholder species names such as *_sp012345678
    df = df[~df["Species"].apply(is_placeholder_species)].copy()

    # Keep only genomes whose fasta files exist
    df["genome_exists"] = df["Genome_ID"].apply(lambda x: genome_file_exists(rep_dir, x))
    df = df[df["genome_exists"]].copy()

    return df


def build_genus_candidates(df: pd.DataFrame) -> dict[tuple[str, str], list[dict]]:
    grouped = defaultdict(list)

    for _, row in df.iterrows():
        key = (row["Domain"], row["Genus"])
        grouped[key].append({
            "Genome_ID": row["Genome_ID"],
            "Species": row["Species"],
            "Genus": row["Genus"],
            "Domain": row["Domain"],
        })

    valid = {}
    for key, records in grouped.items():
        species_set = {r["Species"] for r in records}
        if len(species_set) >= 3:
            valid[key] = records

    return valid


def pick_three_distinct_species(records: list[dict], rng: random.Random) -> list[dict]:
    by_species = defaultdict(list)
    for r in records:
        by_species[r["Species"]].append(r)

    species_list = list(by_species.keys())
    chosen_species = rng.sample(species_list, 3)

    selected = []
    for sp in chosen_species:
        selected.append(rng.choice(by_species[sp]))

    return selected


def choose_four_groups(candidates: dict[tuple[str, str], list[dict]], seed: int) -> list[tuple[tuple[str, str], list[dict]]]:
    rng = random.Random(seed)

    archaea_keys = [k for k in candidates if k[0] == "Archaea"]
    all_keys = list(candidates.keys())

    if len(archaea_keys) < 1:
        raise RuntimeError("No eligible Archaea genus with >=3 non-placeholder species and available genomes was found.")
    if len(all_keys) < 4:
        raise RuntimeError("Fewer than 4 eligible genera available in total after filtering.")

    first_archaea = rng.choice(archaea_keys)
    remaining_keys = [k for k in all_keys if k != first_archaea]
    other_three = rng.sample(remaining_keys, 3)

    chosen_keys = [first_archaea] + other_three
    rng.shuffle(chosen_keys)

    chosen_groups = []
    for key in chosen_keys:
        selected_records = pick_three_distinct_species(candidates[key], rng)
        chosen_groups.append((key, selected_records))

    return chosen_groups


def main() -> None:
    args = parse_args()

    taxonomy_tsv = Path(args.taxonomy_tsv)
    rep_dir = Path(args.rep_dir)

    df = load_taxonomy(taxonomy_tsv, rep_dir)
    candidates = build_genus_candidates(df)
    chosen_groups = choose_four_groups(candidates, args.seed)

    out_rows = []
    print("\nSelected 4 genus-based 3-species groups (non-placeholder species only):\n")
    for idx, ((domain, genus), recs) in enumerate(chosen_groups, start=1):
        print(f"Group {idx}: Domain={domain}, Genus={genus}")
        for r in recs:
            genome_path = rep_dir / f"{r['Genome_ID']}_genomic.fna.gz"
            print(f"  - {r['Genome_ID']}\t{r['Species']}\t{genome_path}")
            out_rows.append({
                "Group": idx,
                "Domain": domain,
                "Genus": genus,
                "Genome_ID": r["Genome_ID"],
                "Species": r["Species"],
                "Genome_File": str(genome_path),
            })
        print()

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.output_tsv, sep="\t", index=False)

    print(f"Saved selection to: {args.output_tsv}")
    print(f"Eligible genus count: {len(candidates)}")
    print(f"Eligible Archaea genus count: {sum(1 for k in candidates if k[0] == 'Archaea')}")


if __name__ == "__main__":
    main()
