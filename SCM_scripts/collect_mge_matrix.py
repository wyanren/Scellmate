#!/usr/bin/env python3
from pathlib import Path
import csv
import argparse

BASE = Path(".").resolve()

def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect idxstats and breadth files into sample-by-MGE matrices"
    )
    parser.add_argument(
        "--prefix",
        default="",
        help='Comma-separated prefixes to keep, e.g. "scDNA_,round_". '
             'Only reference IDs starting with one of these prefixes will be retained.'
    )
    parser.add_argument(
        "--mge-fasta",
        default="",
        help='Comma-separated fasta files containing valid MGE IDs in header lines, '
             'e.g. "../MGE_db/rep-plasmid.fna,../MGE_db/rep-phage.fna"'
    )
    return parser.parse_args()

def get_prefixes(prefix_arg):
    return [x.strip() for x in prefix_arg.split(",") if x.strip()]

def get_mge_fastas(fasta_arg):
    return [x.strip() for x in fasta_arg.split(",") if x.strip()]

def load_mge_ids(fasta_files):
    mge_ids = set()
    for fasta in fasta_files:
        with open(fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    mge_ids.add(line[1:].strip().split()[0])
    return mge_ids

def keep_mge(mge, prefixes, valid_mge_ids):
    if valid_mge_ids is not None:
        return mge in valid_mge_ids
    if not prefixes:
        return True
    return any(mge.startswith(p) for p in prefixes)

def find_files(kind):
    files = []
    for species_dir in BASE.iterdir():
        if not species_dir.is_dir():
            continue
        mapping_dir = species_dir / "mapping"
        if not mapping_dir.exists():
            continue

        if kind == "idxstats":
            candidates = mapping_dir.glob("*.idxstats")
        elif kind == "breadth":
            candidates = mapping_dir.glob("*.breadth")
        else:
            raise ValueError(f"Unknown kind: {kind}")

        for f in candidates:
            name = f.name
            if name.endswith(".raw.idxstats"):
                continue
            if name.endswith(".raw.depth"):
                continue
            if name.endswith(".raw_breadth"):
                continue
            files.append(f)

    return sorted(files)

def parse_idxstats_file(fp, prefixes, valid_mge_ids):
    result = {}
    with open(fp, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            mge = parts[0]
            if mge == "*":
                continue
            if not keep_mge(mge, prefixes, valid_mge_ids):
                continue

            try:
                mapped_reads = float(parts[2])
            except ValueError:
                continue

            result[mge] = mapped_reads
    return result

def parse_breadth_file_absolute(fp, prefixes, valid_mge_ids):
    result = {}
    with open(fp, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            mge = parts[0]
            if not keep_mge(mge, prefixes, valid_mge_ids):
                continue

            try:
                value = float(parts[1])
            except ValueError:
                continue

            result[mge] = value
    return result

def parse_breadth_file_relative(fp, prefixes, valid_mge_ids):
    result = {}
    with open(fp, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            mge = parts[0]
            if not keep_mge(mge, prefixes, valid_mge_ids):
                continue

            try:
                value = float(parts[3])
            except ValueError:
                continue

            result[mge] = value
    return result

def sample_name_from_file(fp, kind):
    name = fp.name
    if kind == "idxstats" and name.endswith(".idxstats"):
        return name[:-9]
    if kind == "breadth" and name.endswith(".breadth"):
        return name[:-8]
    return fp.stem

def build_matrix(files, kind, parser, outname, prefixes, valid_mge_ids):
    all_mges = set()
    per_sample = {}

    for fp in files:
        sample = sample_name_from_file(fp, kind)
        values = parser(fp, prefixes, valid_mge_ids)
        per_sample[sample] = values
        all_mges.update(values.keys())

    all_mges = sorted(all_mges)
    samples = sorted(per_sample.keys())

    with open(outname, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Sample"] + all_mges)
        for sample in samples:
            row = [sample]
            values = per_sample[sample]
            for mge in all_mges:
                row.append(values.get(mge, 0))
            writer.writerow(row)

    print(f"[OK] wrote {outname}")
    print(f"  samples: {len(samples)}")
    print(f"  MGEs:    {len(all_mges)}")

def main():
    args = parse_args()
    prefixes = get_prefixes(args.prefix)
    mge_fastas = get_mge_fastas(args.mge_fasta)
    valid_mge_ids = load_mge_ids(mge_fastas) if mge_fastas else None

    print(f"[INFO] Prefix filter: {prefixes if prefixes else 'None'}")
    print(f"[INFO] MGE fasta files: {mge_fastas if mge_fastas else 'None'}")
    print(f"[INFO] Loaded MGE IDs: {len(valid_mge_ids) if valid_mge_ids is not None else 'None'}")

    idxstats_files = find_files("idxstats")
    breadth_files = find_files("breadth")

    print(f"[INFO] idxstats files found: {len(idxstats_files)}")
    print(f"[INFO] breadth files found: {len(breadth_files)}")

    build_matrix(
        idxstats_files,
        "idxstats",
        parse_idxstats_file,
        "all_samples.idxstats_mapped_reads.matrix.tsv",
        prefixes,
        valid_mge_ids
    )

    build_matrix(
        breadth_files,
        "breadth",
        parse_breadth_file_absolute,
        "all_samples.breadth_absolute.matrix.tsv",
        prefixes,
        valid_mge_ids
    )

    build_matrix(
        breadth_files,
        "breadth",
        parse_breadth_file_relative,
        "all_samples.breadth_relative.matrix.tsv",
        prefixes,
        valid_mge_ids
    )

if __name__ == "__main__":
    main()
