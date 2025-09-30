#!/usr/bin/env python3
"""
Compute subsampling fractions for multiple target read-pair depths.

Inputs:
  - A CSV/TSV file with at least two columns:
      * sample_id          (string)
      * total_readpairs    (integer; total read pairs for that sample)
    Example:
      sample_id,total_readpairs
      S1,42500000
      S2,120000000
      S3,80000000

Usage examples:

  # 1) Auto-generate 10 target depths in millions up to the maximum across samples.
  python subsampling_factors.py samples.tsv --auto 10 --minM 1 --roundM 1 --out all_targets.tsv

  # 2) Provide explicit target depths (in millions):
  python subsampling_factors.py samples.tsv --targets 1,5,10,15,20,25,30,40,50,70 --out all_targets.tsv
  # 3) Also write per-target TSV files, one per target, ready for Nextflow:
  python subsampling_factors.py samples.tsv --targets 1,5,10,15,20,25,30,40,50,70 \
      --per-target-dir factors_by_target/ --out all_targets.tsv

Output (long-format TSV):
  Columns:
    sample_id
    total_readpairs
    target_readpairs
    target_M
    fraction
    estimated_pairs_kept

  Optionally, per-target TSVs in --per-target-dir with columns:
    sample_id, fraction, target_readpairs
"""

import argparse
import math
import os
import sys
from typing import List, Tuple

import pandas as pd


def parse_targets_millions(s: str) -> List[float]:
    vals = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            vals.append(float(tok))
        except ValueError:
            raise SystemExit(f"Could not parse target value '{tok}' as a number.")
    if not vals:
        raise SystemExit("No valid targets parsed.")
    return vals


def auto_targets_millions(k: int, minM: float, maxM: float, roundM: float = 1.0) -> List[float]:
    """
    Create k targets from minM to maxM (inclusive), evenly spaced,
    and rounded to nearest 'roundM' (e.g., 1.0 => nearest 1M).
    Ensures uniqueness and sorted order.
    """
    if k < 2:
        return [round_to(minM, roundM)]
    if maxM <= 0:
        return [round_to(minM, roundM)]
    if minM > maxM:
        minM, maxM = maxM, minM  # swap

    raw = [minM + i * (maxM - minM) / (k - 1) for i in range(k)]
    rounded = [round_to(x, roundM) for x in raw]
    # ensure endpoints are exactly minM and maxM (rounded)
    rounded[0] = round_to(minM, roundM)
    rounded[-1] = round_to(maxM, roundM)

    # deduplicate while preserving order
    uniq = []
    seen = set()
    for x in rounded:
        key = f"{x:.6f}"
        if key not in seen:
            seen.add(key)
            uniq.append(x)

    # if uniqueness collapsed the list too much (e.g., narrow range), fill by stepping by roundM
    if len(uniq) < k:
        step = max(roundM, (maxM - minM) / max(1, (k - 1)))
        uniq = []
        cur = minM
        while cur <= maxM + 1e-9:
            uniq.append(round_to(cur, roundM))
            cur += step
        # dedupe and trim/extend to k length
        uniq2 = []
        seen = set()
        for x in uniq:
            key = f"{x:.6f}"
            if key not in seen:
                seen.add(key)
                uniq2.append(x)
        # ensure min and max present
        if round_to(minM, roundM) not in uniq2:
            uniq2.insert(0, round_to(minM, roundM))
        if round_to(maxM, roundM) not in uniq2:
            uniq2.append(round_to(maxM, roundM))
        uniq = sorted(set(uniq2))
        # If still too many, downsample evenly
        if len(uniq) > k:
            idxs = [round(i * (len(uniq) - 1) / (k - 1)) for i in range(k)]
            uniq = [uniq[i] for i in idxs]
    return sorted(uniq)


def round_to(x: float, base: float) -> float:
    if base <= 0:
        return x
    return round(x / base) * base


def main():
    ap = argparse.ArgumentParser(description="Compute subsampling fractions for multiple target depths.")
    ap.add_argument("input", help="Input CSV/TSV with columns: sample_id,total_readpairs")
    ap.add_argument("--sep", default=None, help="Field separator (default: auto-detect). Use ',' or '\\t' etc.")
    ap.add_argument("--id-col", default="sample_id", help="Column name for sample IDs (default: sample_id)")
    ap.add_argument("--count-col", default="total_readpairs", help="Column name for total read pairs (default: total_readpairs)")

    tgt = ap.add_mutually_exclusive_group()
    tgt.add_argument("--targets", default=None,
                     help="Comma-separated list of target depths in MILLIONS, e.g. '1,5,10,15,20,25,30,40,50,70'")
    tgt.add_argument("--auto", type=int, default=10,
                     help="Auto-generate this many targets (default: 10) spanning minM..maxM in MILLIONS")

    ap.add_argument("--minM", type=float, default=1.0, help="Minimum target (in millions) if using --auto (default: 1.0)")
    ap.add_argument("--roundM", type=float, default=1.0, help="Round targets to nearest X millions (default: 1.0)")
    ap.add_argument("--out", default=None, help="Output TSV path (default: stdout)")
    ap.add_argument("--per-target-dir", default=None,
                    help="If set, also write per-target TSV files into this directory (one file per target).")

    args = ap.parse_args()

    # Read table
    try:
        df = pd.read_csv(args.input, sep=args.sep, engine="python")
    except Exception as e:
        raise SystemExit(f"Failed to read input file: {e}")

    # Normalize column names (case-insensitive match)
    cols_map = {c.lower(): c for c in df.columns}
    if args.id_col.lower() not in cols_map or args.count_col.lower() not in cols_map:
        raise SystemExit(
            f"Input must have columns '{args.id_col}' and '{args.count_col}'. Found columns: {list(df.columns)}"
        )

    id_col = cols_map[args.id_col.lower()]
    count_col = cols_map[args.count_col.lower()]

    # Validate counts
    if df[count_col].isnull().any():
        raise SystemExit(f"Column '{count_col}' contains missing values.")
    try:
        df[count_col] = pd.to_numeric(df[count_col], errors="raise")
    except Exception:
        raise SystemExit(f"Column '{count_col}' must be numeric (total read pairs).")
    if (df[count_col] < 0).any():
        raise SystemExit(f"Column '{count_col}' contains negative values.")
    # Ensure integer-like
    df[count_col] = df[count_col].round().astype(int)

    if df.empty:
        raise SystemExit("Input is empty.")

    # Determine max and target list (in millions)
    max_pairs = int(df[count_col].max())
    maxM = max_pairs / 1e6

    if args.targets is not None:
        targetsM = sorted(set(parse_targets_millions(args.targets)))
    else:
        targetsM = auto_targets_millions(k=args.auto, minM=args.minM, maxM=maxM, roundM=args.roundM)

    # Sanity: filter out <=0 targets
    targetsM = [t for t in targetsM if t > 0]
    if not targetsM:
        raise SystemExit("No valid target depths after processing.")

    # Compute
    out_rows = []
    for _, row in df.iterrows():
        sid = str(row[id_col])
        total = int(row[count_col])
        for tM in targetsM:
            target_pairs = int(round(tM * 1_000_000))
            if total <= 0:
                frac = 0.0
                kept = 0
            else:
                frac = min(1.0, target_pairs / total)
                kept = int(math.floor(total * frac))
            # format fraction to reasonable precision for tools like seqtk
            frac_fmt = float(f"{frac:.8f}")
            out_rows.append(
                {
                    "sample_id": sid,
                    "total_readpairs": total,
                    "target_readpairs": target_pairs,
                    "target_M": tM,
                    "fraction": frac_fmt,
                    "estimated_pairs_kept": kept,
                }
            )

    out_df = pd.DataFrame(out_rows).sort_values(["target_readpairs", "sample_id"]).reset_index(drop=True)

    # Write main output
    if args.out:
        out_df.to_csv(args.out, sep="\t", index=False)
    else:
        out_df.to_csv(sys.stdout, sep="\t", index=False)

    # Optional: per-target files for easy use in Nextflow
    if args.per_target_dir:
        os.makedirs(args.per_target_dir, exist_ok=True)
        for tM in sorted(set(out_df["target_M"])):
            tsub = out_df[out_df["target_M"] == tM].copy()
            tsub = tsub[["sample_id", "fraction", "target_readpairs"]]
            fname = os.path.join(args.per_target_dir, f"subsample_factors_T{int(round(tM))}M.tsv")
            tsub.to_csv(fname, sep="\t", index=False)


if __name__ == "__main__":
    main()
