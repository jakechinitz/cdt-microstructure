#!/usr/bin/env python3
"""Is EPRL's high blob (~5.0) a REAL sharper-de-Sitter feature, or an artifact?
And does it square with the d_s(sigma) flow vs Regge?

blob_score = max(V3)/mean(V3) is a CRUDE peakedness ratio, NOT a de Sitter test:
a spike (one fat slice, rest stalk) inflates it identically to a smooth cos^3 bump.
The de Sitter content is elsewhere:
  * cos3_relerr  LOW  => profile actually fits A*cos^3 (de Sitter-like)
                 HIGH => high max/mean but NOT cos^3 (spike artifact)
  * active_slices ~several => blob spans slices (de Sitter); ~1 => a spike
  * d_s(sigma) flow => the DIRECT manifold-vs-not channel
Verdict: high blob + LOW cos3_relerr + several active = genuinely sharper de Sitter.
         high blob + HIGH cos3_relerr + ~1 active   = spike, no prediction.

Two more cautions baked in:
  - EPRL may not be N4-equilibrated (its blob/d_s are then transient). This script
    flags whether each run's N4 is still trending over the tail.
  - The d_s comparison is done at MATCHED N4 (nearest curve in each sidecar to a
    target), since d_s reads the N4-sized dual graph. Caveat printed: a run that
    reached the target N4 while still climbing is less thermalized there than one
    equilibrated at it.

Usage (from v6_vm_bundle/):
  python blob_ds_diagnose.py                 # default target N4 = 85000 (overlap top)
  python blob_ds_diagnose.py --target-n4 80000 --tail 8
"""
from __future__ import annotations
import argparse
import json
import re
import numpy as np

# log columns: sweep N4 N41 d_H blob active cos3err valid links wall_s
ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+"
                 r"([\d.]+|--)\s+(\w+)")


def log_tail(path, tail):
    rows = []
    try:
        for ln in open(path):
            m = ROW.match(ln)
            if m:
                cos = m[7]
                rows.append(dict(sweep=int(m[1]), N4=int(m[2]), N41=int(m[3]),
                                 dH=float(m[4]), blob=float(m[5]),
                                 active=int(m[6]),
                                 cos3=(float(cos) if cos != "--" else None),
                                 valid=m[8]))
    except FileNotFoundError:
        return None
    return rows[-tail:] if rows else []


def n4_trending(rows):
    n4 = np.array([r["N4"] for r in rows], float)
    sw = np.array([r["sweep"] for r in rows], float)
    if len(n4) < 3:
        return True
    drift = abs(np.polyfit(sw - sw.mean(), n4, 1)[0]) * (sw[-1] - sw[0])
    return drift > 0.02 * n4.mean()


def load_curves(path):
    out = []
    try:
        for ln in open(path):
            ln = ln.strip()
            if ln:
                out.append(json.loads(ln))
    except FileNotFoundError:
        return None
    return out


def nearest_curve(curves, target_n4):
    return min(curves, key=lambda c: abs(c["N4"] - target_n4)) if curves else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-n4", type=int, default=85000)
    ap.add_argument("--tail", type=int, default=8)
    args = ap.parse_args()

    print("=" * 70)
    print("PART 1 -- is EPRL's high blob real de Sitter or a spike? (blob is max/mean,")
    print("         NOT a de Sitter test; cos3_relerr + active_slices decide)")
    print("=" * 70)
    print(f"{'run':6} {'N4':>8} {'blob':>6} {'cos3err':>8} {'active':>7} {'N4-trending?':>13}")
    diag = {}
    for k, path in (("floor", "logs/ladder_floor.log"),
                    ("eprl", "logs/ladder_eprl.log"),
                    ("regge", "logs/ladder_regge.log")):
        rows = log_tail(path, args.tail)
        if not rows:
            print(f"{k:6}  (no data)")
            continue
        blob = np.median([r["blob"] for r in rows])
        cos = [r["cos3"] for r in rows if r["cos3"] is not None]
        cos_m = np.median(cos) if cos else float("nan")
        act = np.median([r["active"] for r in rows])
        tr = n4_trending(rows)
        diag[k] = dict(blob=blob, cos3=cos_m, active=act, trending=tr,
                       N4=rows[-1]["N4"])
        print(f"{k:6} {rows[-1]['N4']:>8} {blob:>6.2f} {cos_m:>8.3f} {act:>7.0f} "
              f"{'YES (transient!)' if tr else 'no':>13}")

    if "eprl" in diag:
        e = diag["eprl"]
        print("\n  EPRL blob verdict:")
        if e["cos3"] != e["cos3"]:  # nan
            print("    cos3_relerr unavailable -- can't classify.")
        elif e["cos3"] < 0.5 and e["active"] >= 3:
            print(f"    blob={e['blob']:.1f} with LOW cos3err={e['cos3']:.3f} and "
                  f"{e['active']:.0f} active slices => looks like GENUINE de Sitter,")
            print("    not a spike. (Still must clear the equilibration caveat below.)")
        else:
            print(f"    blob={e['blob']:.1f} but cos3err={e['cos3']:.3f} "
                  f"(active={e['active']:.0f}) => NOT a clean cos^3. The high blob is")
            print("    likely a SPIKE / profile artifact, not 'sharper de Sitter'. Do")
            print("    NOT publish it as a de Sitter prediction.")
        if e["trending"]:
            print("    !! EPRL N4 STILL CLIMBING -> blob is a TRANSIENT measurement;")
            print("       it may change as N4 equilibrates. Not bankable yet.")

    print("\n" + "=" * 70)
    print(f"PART 2 -- d_s(sigma) flow, EPRL vs Regge at matched N4~{args.target_n4}")
    print("         (the DIRECT manifold channel; resolves the blob-vs-d_s tension)")
    print("=" * 70)
    ce = load_curves("results/ds_eprl.jsonl")
    cr = load_curves("results/ds_regge.jsonl")
    if not ce or not cr:
        print("  (missing d_s sidecar(s))")
        return
    qe, qr = nearest_curve(ce, args.target_n4), nearest_curve(cr, args.target_n4)
    print(f"  EPRL  curve at N4={qe['N4']}  (target {args.target_n4})")
    print(f"  Regge curve at N4={qr['N4']}")
    if abs(qe["N4"] - qr["N4"]) > 0.1 * args.target_n4:
        print(f"  !! nearest curves differ in N4 by {abs(qe['N4']-qr['N4'])} "
              f"(>10%) -- comparison is N4-confounded, treat as indicative only.")
    se, de = np.array(qe["sigma"]), np.array(qe["d_s"])
    sr, dr = np.array(qr["sigma"]), np.array(qr["d_s"])
    # print the two flows at a few sigma bands
    print(f"\n  {'sigma':>7} {'d_s(EPRL)':>10} {'d_s(Regge)':>11} {'diff':>7}")
    idxs = np.linspace(0, len(se) - 1, 8).astype(int)
    for i in idxs:
        # match regge sigma index (same logspace grid)
        j = min(range(len(sr)), key=lambda x: abs(sr[x] - se[i]))
        print(f"  {se[i]:>7.2f} {de[i]:>10.2f} {dr[j]:>11.2f} {de[i]-dr[j]:>+7.2f}")
    pe_lo = np.median(de[se < np.percentile(se, 25)])
    pe_hi = np.median(de[se > np.percentile(se, 60)])
    pr_lo = np.median(dr[sr < np.percentile(sr, 25)])
    pr_hi = np.median(dr[sr > np.percentile(sr, 60)])
    print(f"\n  EPRL  flow: small-sigma {pe_lo:.2f} -> large-sigma {pe_hi:.2f}  "
          f"(rise {pe_hi-pe_lo:+.2f})")
    print(f"  Regge flow: small-sigma {pr_lo:.2f} -> large-sigma {pr_hi:.2f}  "
          f"(rise {pr_hi-pr_lo:+.2f})")
    print("\n  READ: if EPRL's large-sigma d_s and rise ~ Regge's, then d_s says")
    print("  'same as Regge' while blob says 'sharper' -> the blob is the OUTLIER and")
    print("  PART 1's cos3err/active is the tiebreaker (spike vs real). If EPRL's d_s")
    print("  rise EXCEEDS Regge's too, blob and d_s agree on 'sharper' -- a real (if")
    print("  still equilibration-caveated) signal worth taking seriously.")


if __name__ == "__main__":
    main()
