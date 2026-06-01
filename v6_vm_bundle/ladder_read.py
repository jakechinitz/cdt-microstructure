#!/usr/bin/env python3
"""Ladder read: GATE 0 (equilibration overlap) -> separation -> eprl position.

Reads the three d_s(sigma) sidecars (ds_floor / ds_eprl / ds_regge .jsonl), each a
JSONL of {sweep, N4, N41, d_H, sigma[], d_s[]} per measurement, and applies the
pre-registered read ORDER. It refuses to advance a gate until the prior one passes.

GATE 0 -- equilibration overlap (the load-bearing precondition):
  The three actions pin N41 but settle at DIFFERENT N4 (no_action high via N32
  entropy; Regge low via its k4 term). d_s reads the N4-sized dual graph, so a
  comparison is only valid at a common N4 where BOTH runs are near-equilibrated.
  This gate finds, for each pair, an N4 window where both runs are equilibrated
  (d_s stable across sweeps passing through that N4). If no such window exists,
  the comparison is NOT readable -- the honest finding is "N41-pinned CDT cannot
  produce volume-and-equilibration-matched rungs for actions with divergent volume
  preferences" (an instrument-limit result). Saturating-high is just as fatal as
  runaway if high does not overlap low.

GATE 1 -- separation: at a Gate-0-valid common N4, do floor and Regge d_s(sigma)
  FLOW SHAPES resolve apart? If not -> "can't resolve gravity from expander at
  feasible volume", and eprl is unreadable.

GATE 2 -- eprl position: only if Gates 0+1 pass. Where does eprl's d_s(sigma) sit
  between floor and Regge at the common N4? Near-floor = null (pre-registered
  likely); near-Regge = surprise/signal. Claim: structural compatibility on
  CDT-proposed geometry, NOT generation.

Usage (from v6_vm_bundle/, after the runs have grown):
  python ladder_read.py
  python ladder_read.py --n4-tol 0.08 --equil-cv 0.04
"""
from __future__ import annotations
import argparse
import json
import numpy as np


def load(path):
    rows = []
    try:
        for ln in open(path):
            ln = ln.strip()
            if ln:
                rows.append(json.loads(ln))
    except FileNotFoundError:
        return None
    return rows


def equilibrated_windows(rows, equil_cv=0.04, min_pts=3):
    """Return list of (N4, d_s_curve, sweep) for measurements that are
    'equilibrated' = the d_H is stable (low coefficient of variation) over a
    trailing window of measurements around it. Crude but data-driven."""
    if not rows or len(rows) < min_pts:
        return []
    dH = np.array([r["d_H"] for r in rows])
    out = []
    for i in range(min_pts - 1, len(rows)):
        seg = dH[i - min_pts + 1:i + 1]
        cv = seg.std() / max(abs(seg.mean()), 1e-9)
        if cv < equil_cv:
            out.append(rows[i])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--floor", default="results/ds_floor.jsonl")
    ap.add_argument("--eprl", default="results/ds_eprl.jsonl")
    ap.add_argument("--regge", default="results/ds_regge.jsonl")
    ap.add_argument("--n4-tol", type=float, default=0.08,
                    help="two runs are at 'common N4' if their N4 differ < this frac")
    ap.add_argument("--equil-cv", type=float, default=0.04,
                    help="d_H coeff-of-variation below which a window counts as equilibrated")
    args = ap.parse_args()

    runs = {"floor": load(args.floor), "eprl": load(args.eprl), "regge": load(args.regge)}
    for k, v in runs.items():
        if not v:
            print(f"[wait] {k}: no data yet ({args.floor if k=='floor' else args.eprl if k=='eprl' else args.regge})")
    if not all(runs.values()):
        print("Not all sidecars have data; let the runs grow.")
        return

    print("=== run summary ===")
    for k, v in runs.items():
        n4 = [r["N4"] for r in v]
        n41 = [r["N41"] for r in v]
        print(f"  {k:6}: {len(v)} curves, N4 {min(n4)}..{max(n4)}, "
              f"N41~{int(np.median(n41))}, N4/N41~{np.median(n4)/np.median(n41):.1f}")

    eq = {k: equilibrated_windows(v, args.equil_cv) for k, v in runs.items()}
    print("\n=== GATE 0: equilibration overlap ===")
    for k in runs:
        if eq[k]:
            n4s = [r["N4"] for r in eq[k]]
            print(f"  {k:6}: equilibrated N4 range {min(n4s)}..{max(n4s)} ({len(eq[k])} pts)")
        else:
            print(f"  {k:6}: NO equilibrated window yet (d_H still drifting)")

    def overlap(a, b):
        if not eq[a] or not eq[b]:
            return None
        na = [r["N4"] for r in eq[a]]; nb = [r["N4"] for r in eq[b]]
        lo, hi = max(min(na), min(nb)), min(max(na), max(nb))
        return (lo, hi) if lo <= hi else None

    pairs = [("floor", "regge"), ("floor", "eprl"), ("eprl", "regge")]
    gate0 = {}
    for a, b in pairs:
        ov = overlap(a, b)
        gate0[(a, b)] = ov
        if ov:
            print(f"  OVERLAP {a}-{b}: common equilibrated N4 in [{ov[0]}, {ov[1]}] -> READABLE")
        else:
            print(f"  OVERLAP {a}-{b}: NONE -> NOT readable at matched equilibrated volume")

    if not gate0[("floor", "regge")]:
        print("\n[GATE 0 FAIL on floor-regge] The ladder's defining pair has no common")
        print("  equilibrated N4. HONEST FINDING: N41-pinned CDT cannot produce volume-")
        print("  and-equilibration-matched rungs for actions with divergent volume")
        print("  preferences. Instrument-limit result. d_s comparison NOT readable;")
        print("  do NOT read eprl against either endpoint.")
        return

    print("\n=== GATE 1: floor-vs-regge separation (at common equilibrated N4) ===")
    print("  [run again once both have equilibrated curves in the overlap window;]")
    print("  [compare d_s(sigma) flow SHAPE there -- does Regge bend up at large sigma]")
    print("  [while floor stays flat? If not -> resolution-ceiling finding.]")
    print("\n=== GATE 2: eprl position (only if Gates 0+1 pass) ===")
    print("  [near-floor = null (pre-registered likely); near-Regge = signal.]")
    print("  [claim: structural compatibility on CDT-proposed geometry, NOT generation.]")


if __name__ == "__main__":
    main()
