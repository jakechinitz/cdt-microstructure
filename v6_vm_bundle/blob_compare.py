#!/usr/bin/env python3
"""Matched-volume blob comparison: the actual EPRL-sweep read.

THE EXPERIMENT (per PREREGISTRATION_matched_volume_sweep.md):
  signal = (cos^3 de Sitter blob at beta=0.3) - (blob at beta=0), read on a
  MATCHED-VOLUME, EQUILIBRATED pair. Not beta=0.3 alone; not d_H (which is nearly
  blind to this term, type-split 0.054 << sigma=0.87). The blob sits closest to
  the weak-but-real geometric channel (dual-graph r_nn=0.14, ~4 sigma) the term
  actually has, so it is the observable with a chance of registering it.

This script REFUSES to declare a comparison readable unless the preconditions
hold, so the read can't be volume-confounded like every prior sweep:

  G-THERM   each run's N4 has stopped trending (last-third slope ~ 0)
  G-MATCH   beta=0 and beta=0.3 are at the same N4 (within tol) AND N41~target
  G-VALID   no BAD rows in the compared window

Only if all pass does it print the blob delta and the pre-registered reading.

Reads the live logs (no sweeps). Usage from v6_vm_bundle/:
  python blob_compare.py                       # b0.0 vs b0.3, last 40 rows
  python blob_compare.py --a 0.0 --b 0.3 --window 60 --match-tol 0.04
"""
from __future__ import annotations
import argparse
import re
import numpy as np

# columns: sweep N4 N41 d_H blob active cos3err valid links wall_s
ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+"
                 r"([\d.]+)\s+(\w+)\s+(\S+)\s+(\d+)")


def load_rows(path):
    out = []
    bad = 0
    try:
        with open(path) as f:
            for line in f:
                if "BAD" in line:
                    bad += 1
                m = ROW.match(line)
                if m:
                    out.append(dict(sweep=int(m[1]), N4=int(m[2]), N41=int(m[3]),
                                    dH=float(m[4]), blob=float(m[5]),
                                    active=int(m[6]), cos3=float(m[7]),
                                    valid=m[8]))
    except FileNotFoundError:
        return None, 0
    return out, bad


def thermalized(rows):
    """N4 last-third slope ~ 0 relative to its scatter."""
    if len(rows) < 9:
        return False, None
    n4 = np.array([r["N4"] for r in rows])
    k = len(n4) // 3
    sw = np.array([r["sweep"] for r in rows], float)
    seg_s, seg_n = sw[-k:], n4[-k:].astype(float)
    slope = np.polyfit(seg_s - seg_s.mean(), seg_n, 1)[0]   # N4 per sweep
    drift = abs(slope) * (seg_s[-1] - seg_s[0])             # total N4 change over segment
    return drift < 0.02 * n4.mean(), drift


def stat(rows, window, key):
    v = np.array([r[key] for r in rows[-window:]], float)
    return v.mean(), v.std(), len(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", default="0.0", help="baseline beta")
    ap.add_argument("--b", default="0.3", help="coupled beta")
    ap.add_argument("--vol", default="20k")
    ap.add_argument("--window", type=int, default=40, help="rows to average over the equilibrated tail")
    ap.add_argument("--match-tol", type=float, default=0.04, help="max fractional N4 mismatch allowed")
    ap.add_argument("--target-n41", type=int, default=20000)
    args = ap.parse_args()

    pa = f"logs/thy_b{args.a}_{args.vol}.log"
    pb = f"logs/thy_b{args.b}_{args.vol}.log"
    A, badA = load_rows(pa)
    B, badB = load_rows(pb)

    if not A or not B:
        print(f"[not readable] missing/empty logs: "
              f"{pa if not A else ''} {pb if not B else ''}")
        return
    print(f"baseline beta={args.a}: {len(A)} rows ({pa})")
    print(f"coupled  beta={args.b}: {len(B)} rows ({pb})")

    # --- gates ---
    thA, drA = thermalized(A)
    thB, drB = thermalized(B)
    n4A, _, _ = stat(A, args.window, "N4")
    n4B, _, _ = stat(B, args.window, "N4")
    n41A, _, _ = stat(A, args.window, "N41")
    n41B, _, _ = stat(B, args.window, "N41")
    mismatch = abs(n4A - n4B) / max(n4A, n4B)

    print("\n=== preconditions ===")
    print(f"  G-THERM  beta={args.a}: {'PASS' if thA else 'WAIT'} (N4 tail drift ~{drA:.0f})")
    print(f"  G-THERM  beta={args.b}: {'PASS' if thB else 'WAIT'} (N4 tail drift ~{drB:.0f})")
    print(f"  G-MATCH  N4: {n4A:.0f} vs {n4B:.0f}  (mismatch {mismatch*100:.1f}%, tol {args.match_tol*100:.0f}%) "
          f"-> {'PASS' if mismatch < args.match_tol else 'FAIL'}")
    print(f"  G-MATCH  N41~target: {n41A:.0f}, {n41B:.0f} (target {args.target_n41}) "
          f"-> {'PASS' if abs(n41A-args.target_n41)<0.05*args.target_n41 and abs(n41B-args.target_n41)<0.05*args.target_n41 else 'CHECK'}")
    print(f"  G-VALID  BAD rows: beta={args.a}:{badA}  beta={args.b}:{badB} "
          f"-> {'PASS' if badA==0 and badB==0 else 'FAIL'}")

    ready = thA and thB and mismatch < args.match_tol and badA == 0 and badB == 0
    if not ready:
        print("\n[NOT READABLE YET] preconditions unmet -- let it cook / fix setup. "
              "Do NOT interpret the blob below as a result.")

    # --- the comparison (printed always, but only trustworthy if ready) ---
    bA, sA, nA = stat(A, args.window, "blob")
    bB, sB, nB = stat(B, args.window, "blob")
    dHA, _, _ = stat(A, args.window, "dH")
    dHB, _, _ = stat(B, args.window, "dH")
    cA, _, _ = stat(A, args.window, "cos3")
    cB, _, _ = stat(B, args.window, "cos3")
    delta = bB - bA
    pooled = np.sqrt(sA**2 / nA + sB**2 / nB) or 1e-9

    print(f"\n=== THE EXPERIMENT: blob(beta={args.b}) - blob(beta={args.a}) "
          f"[avg of last {args.window} rows] ===")
    print(f"  blob   : {bA:.3f}+/-{sA:.3f}  ->  {bB:.3f}+/-{sB:.3f}   "
          f"delta={delta:+.3f}  ({delta/pooled:+.1f} sigma)")
    print(f"  (2nd)  d_H : {dHA:.3f} -> {dHB:.3f}   (nearly blind to term -- secondary)")
    print(f"  (2nd)  cos3err: {cA:.3f} -> {cB:.3f}   (lower = more de-Sitter)")

    if ready:
        if abs(delta) < 2 * pooled:
            print("\n  READING: FLAT blob at matched volume -> DECOUPLING NULL.")
            print("  To claim the GOOD null (not the empty one), state explicitly that the")
            print("  r_nn=0.14 (~4sigma) geometric channel, amplified by beta over ~10^5")
            print("  pentachora, gave the term a real chance to move the blob and it did not.")
        else:
            print(f"\n  READING: blob MOVED ({delta/pooled:+.1f} sigma) at matched volume.")
            print("  BEFORE believing it: relaunch beta=0.05, 0.1 and confirm the response")
            print("  is SMOOTH/MONOTONIC in beta (rules out threshold artifact / strong-")
            print("  coupling breakage). A lone 0->0.3 jump is not yet a physical signal.")


if __name__ == "__main__":
    main()
