#!/usr/bin/env python3
"""
run_spades_strict.py — Parallel SPAdes runner (fixed parameters)

Features
--------
1. **Sample count announcement** before execution starts.
2. **Live progress bar** based on `tqdm`. Falls back to simple counter if `tqdm` is not available.
3. **Silent SPAdes execution**: all STDOUT/STDERR from every SPAdes run are redirected to a single `spades.log` inside the output directory (`-o`).

Example
-------
python run_spades_strict.py \
    -i /path/to/01_trim_SAGs \ # the diretory of paired reads
    -o /path/to/04_SAG_assembly/spades_output \
    --suffix_R1 _R1_paired.fastq \
    --threads 40 --per-job 4
"""
import argparse
import concurrent.futures as cf
import pathlib
import shutil
import subprocess
import sys
from datetime import datetime

try:
    from tqdm import tqdm  # type: ignore
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# ──────────────────────────── CLI ────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--indir", required=True, type=pathlib.Path,
                   help="Directory containing R1 fastq files")
    p.add_argument("-o", "--outdir", required=True, type=pathlib.Path,
                   help="Root output directory for SPAdes results")
    p.add_argument("--suffix_R1", default="_R1_paired.fastq",
                   help="Suffix of R1 files used to derive sample names")
    p.add_argument("--threads", type=int, default=40,
                   help="Total available CPU threads")
    p.add_argument("--per-job", type=int, default=4,
                   help="Threads per SPAdes process")
    return p.parse_args()

# ────────── run SPAdes for a single sample ──────────

def run_one(args):
    sample, r1, r2, sample_outdir, per_job, root_outdir, log_path = args
    sample_outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "spades.py",
        "-t", str(per_job),
        "--sc", "--careful",
        "--pe1-s", str(r1),
        "--pe1-s", str(r2),
        "-o", str(sample_outdir)
    ]

    start_time = datetime.now().strftime("%F %T")
    with open(log_path, "a") as logf:
        logf.write(f"\n[{start_time}] START {sample}\n")
        ret = subprocess.run(cmd, stdout=logf, stderr=logf).returncode
        end_time = datetime.now().strftime("%F %T")
        logf.write(f"[{end_time}] END   {sample} (exit {ret})\n")

    if ret == 0:
        src = sample_outdir / "contigs.fasta"
        dst = root_outdir / f"{sample}.fasta"
        try:
            shutil.move(src, dst)
        except FileNotFoundError:
            with open(log_path, "a") as logf:
                logf.write(f"[WARN] {src} not found; skipping move\n")
        # Remove intermediate directory to save space
        shutil.rmtree(sample_outdir, ignore_errors=True)
    return ret

# ──────────────────────────── main ───────────────────────────

def main():
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    r1_files = sorted(a.indir.glob(f"*{a.suffix_R1}"))
    if not r1_files:
        sys.exit(f"[ERR] No files matching *{a.suffix_R1} in {a.indir}")

    tasks = []
    for r1 in r1_files:
        sample = r1.name.replace(a.suffix_R1, "")
        r2 = r1.with_name(r1.name.replace("_R1", "_R2"))
        if not r2.exists():
            print(f"[WARN] Missing R2 for {sample}; skipping")
            continue
        sample_outdir = a.outdir / f"{sample}_scCareful"
        tasks.append((sample, r1, r2, sample_outdir, a.per_job, a.outdir, a.outdir / "spades.log"))

    total = len(tasks)
    if total == 0:
        sys.exit("[ERR] No valid sample pairs found.")

    jobs = max(1, a.threads // a.per_job)
    print(f"[INFO] Total samples: {total} | Workers: {jobs} | per_job: {a.per_job}")

    progress_iter = tqdm(total=total, desc="SPAdes", unit="sample") if TQDM_AVAILABLE else None

    fails = 0
    with cf.ProcessPoolExecutor(max_workers=jobs) as pool:
        for ret in pool.map(run_one, tasks):
            fails += (ret != 0)
            if progress_iter:
                progress_iter.update(1)
            else:
                done = total - fails if fails else total  # approximation
                print(f"[INFO] Finished {done}/{total} samples", end="\r", flush=True)

    if progress_iter:
        progress_iter.close()

    print("\n[DONE] " + ("All succeeded" if fails == 0 else f"{fails} samples failed"))
    sys.exit(fails)

if __name__ == "__main__":
    main()

