#!/usr/bin/env python3
"""Ladder read: Gate-0-INDEPENDENT signals first, THEN gated cross-run comparison.

Reads the three d_s(sigma) sidecars (ds_floor / ds_eprl / ds_regge .jsonl), each a
JSONL of {sweep, N4, N41, d_H, sigma[], d_s[]} per measurement.

IMPORTANT (do not let Gate 0 eat the run): Gate 0 only blocks ONE comparison --
the thermalized-endpoint cross-run d_s at matched volume. Several OTHER signals are
readable regardless of whether Gate 0 passes, and are printed FIRST,
unconditionally:
  SIGNAL A -- N4/N41 ratio ordering: what each action does to volume composition at
    fixed N41 (no_action high via N32 entropy; Regge low via k4). EPRL near Regge =
    GR-like; near no_action = entropic. The single most useful readable result, and
    it needs no matching.
  SIGNAL B -- per-run d_s(sigma) flow SHAPE: does each run's own flow bend up
    (manifold) or stay flat (expander)? Characterizes each curve without cross-run
    matching.
  SIGNAL C -- viability: does eprl_only hold a stable valid manifold / pin N41 at
    all under pure-amplitude weighting?
These reach the COMPATIBILITY claim (j=3 structure is CDT-like vs expander-like),
NOT generation. Modest-but-real, and robust to Gate 0 failing.

THEN the gated cross-run comparison (the crisp result, only if matchable):

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


def equilibrated_windows(rows, equil_cv=0.04, min_pts=4):
    """Return measurements that are equilibrated = N4 has STOPPED TRENDING over a
    trailing window. CRITICAL: test N4, not d_H. d_H can plateau while N4 is still
    climbing monotonically (observed: EPRL d_H~3.25 flat while N4 climbs 93k->99k),
    which would falsely label a still-growing run 'equilibrated'. Since the read
    requires matched VOLUME (N4), equilibration must be N4-stability: the trailing
    N4 slope must be ~0 relative to its own scatter."""
    if not rows or len(rows) < min_pts:
        return []
    out = []
    for i in range(min_pts - 1, len(rows)):
        seg = rows[i - min_pts + 1:i + 1]
        n4 = np.array([r["N4"] for r in seg], float)
        sw = np.array([r["sweep"] for r in seg], float)
        # net drift over the window vs the window's own N4 scale
        slope = np.polyfit(sw - sw.mean(), n4, 1)[0]
        drift = abs(slope) * (sw[-1] - sw[0])
        if drift < equil_cv * n4.mean():
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

    # ----------------------------------------------------------------------
    # READABLE SIGNALS THAT DO NOT DEPEND ON GATE 0 (cross-run matching).
    # These survive an overlap failure -- print them FIRST and unconditionally,
    # so an instrument-limit Gate-0 failure never eats the rest of the run.
    # ----------------------------------------------------------------------
    print("\n=== SIGNAL A: N4/N41 ratio ordering (Gate-0-INDEPENDENT) ===")
    print("  What each action does to volume COMPOSITION at fixed N41. No matching")
    print("  needed. no_action high (entropy maxes N32); Regge low (k4 suppresses N4).")
    ratios = {}
    trending = {}
    for k, v in runs.items():
        tail = v[-5:] if len(v) >= 5 else v
        r = float(np.median([row["N4"] / max(row["N41"], 1) for row in tail]))
        ratios[k] = r
        # is N4 still climbing across the tail? then this ratio is a TRANSIENT.
        n4 = np.array([row["N4"] for row in tail], float)
        sw = np.array([row["sweep"] for row in tail], float)
        drift = abs(np.polyfit(sw - sw.mean(), n4, 1)[0]) * (sw[-1] - sw[0])
        trending[k] = drift > 0.02 * n4.mean()
        flag = "  <-- STILL CLIMBING (transient, NOT an equilibrium value)" if trending[k] else ""
        print(f"  {k:6}: N4/N41 = {r:.2f}{flag}")
    if trending.get("eprl"):
        print("  !! EPRL N4 has NOT equilibrated -- its ratio is mid-flight, climbing")
        print("     through Regge's value toward an unknown ceiling. DO NOT bank")
        print("     'EPRL ~ Regge'; it may keep rising toward the entropic floor.")
    if all(k in ratios for k in ("floor", "eprl", "regge")):
        df, dr = abs(ratios["eprl"] - ratios["floor"]), abs(ratios["eprl"] - ratios["regge"])
        verdict = ("CDT/Regge-like (GR-like volume composition)" if dr < df else
                   "expander/no-action-like (entropic volume composition)" if df < dr else
                   "intermediate")
        print(f"  -> EPRL ratio {ratios['eprl']:.2f} is closer to "
              f"{'Regge' if dr < df else 'no_action' if df < dr else 'neither'} "
              f"=> j=3 structure looks {verdict}")
        print("  (compatibility channel, NOT generation; readable regardless of Gate 0.)")

    print("\n=== SIGNAL B: per-run d_s(sigma) FLOW SHAPE (Gate-0-INDEPENDENT) ===")
    print("  Each run's OWN flow character: does d_s bend up at large sigma (manifold)")
    print("  or stay flat (expander)? No cross-run common-N4 needed.")
    for k, v in runs.items():
        last = v[-1]
        ds = np.array(last["d_s"]); sig = np.array(last["sigma"])
        # short-scale vs large-scale d_s on this run's own latest curve
        lo = float(np.median(ds[sig < np.percentile(sig, 25)]))
        hi = float(np.median(ds[sig > np.percentile(sig, 60)]))
        shape = ("MANIFOLD-like (bends up at large scale)" if hi > lo + 0.3 else
                 "FLAT/expander-like" if abs(hi - lo) <= 0.3 else "falling")
        print(f"  {k:6} (N4={last['N4']}): d_s small-σ~{lo:.2f} large-σ~{hi:.2f} -> {shape}")
    print("  -> compare EPRL's shape to floor vs Regge: does it have the manifold")
    print("     character (like Regge) or the flat character (like no_action)?")

    print("\n=== SIGNAL C: viability (Gate-0-INDEPENDENT) ===")
    for k, v in runs.items():
        n41 = [r["N41"] for r in v]
        held = max(abs(np.array(n41) - np.median(n41))) < 0.1 * np.median(n41)
        print(f"  {k:6}: N41 held={held}, grew to N4={v[-1]['N4']} "
              f"(valid manifold throughout = run didn't abort)")
    print("  -> eprl_only holding a stable valid manifold at all = basic geometric")
    print("     viability of pure-amplitude weighting, independent of any matching.")

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
