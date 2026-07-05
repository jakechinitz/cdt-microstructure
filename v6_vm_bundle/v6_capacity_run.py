#!/usr/bin/env python3
"""STAGE 4 -- capacity propagation: does the theory's conserved strain field
carry a defect's influence to distance, with the right form?

Implements Model P from CAPACITY_CONSERVATION.md (Route B: the conservation
reading, the only one that can produce long-range gravity). The design goal
is to measure the two quantities that are NOT preordained:

  M1  EMERGENT SINK LINEARITY. The sink is NOT inserted: cells absorb
      capacity at a rate proportional to their instantaneous closure
      failure (label collisions), which the defect induces in itself and
      its neighborhood. Pinning defects at four commitment levels
      (1/2/3/6 colliding pairs -- a mass ladder) and measuring the flux
      each actually draws tests the paper's 'maintenance linear in
      commitment' postulate dynamically.
  M2  THE JUNCTION CONSTANT. Far-field deficit amplitude per unit of
      measured sink flux -- the lattice analog of Newton's constant,
      produced by the microphysics rather than chosen.

The 1/r form itself is a theorem given conservation (capacity_lattice_law
demo); here it appears as the T3/T4 GATES, not as a discovery:
  T1  exact conservation: sum(f) constant to machine precision every sweep
      (transport is antisymmetric by construction; absorption is recycled
      uniformly -- the closed-manifold analog of 'returned at infinity').
  T2  vacuum anchor: field initialized and maintained at <f> = 7.4198.
  T3  deficit profile around pins matches the massless graph Green
      function (1/r-like on large slices).
  T4  screening-mass fit consistent with zero; a nonzero mass means
      conservation is leaking (Route-A contamination) and voids the run.

MODEL P (declared forks):
  * field f_x >= 0 on spatial cells, one field per slice (the strain field
    of the paper's static branch is spatial); transport ONLY within slices.
  * per sweep: (i) one label heat-bath pass (beta-consistent, pins frozen);
    (ii) transport f' = f + D * sum_faces (f_y - f_x), D <= 1/4;
    (iii) absorption dA_x = kappa * rate_x * f_x, recycled uniformly
    within the slice the same sweep (steady state possible on a closed
    manifold). In --absorb excess mode the rate is computed from the
    PERSISTENT failure average A_x (EWMA over --persist sweeps):
    rate = max(0, A_x - nbar - margin). Rationale (v1.2): the vacuum's
    instantaneous n_coll is ~92% zero / ~8% one, so max(0, n_coll - nbar)
    rectifies the fluctuations and removes only ~8% of the vacuum
    absorption (rect = 0.071 vs nbar = 0.077 at the theory point) --
    giving the vacuum an intrinsic screening length
    xi = sqrt(D/(kappa*rect)) ~ 7.5 steps BY CONSTRUCTION. The v1.1 T4
    'screening' verdict measured that implementation artifact, not the
    theory. The paper defines mass as PERSISTENT closure failure, so
    maintenance is drawn on the time-averaged failure: the window
    suppresses the rectified vacuum rate ~1/sqrt(W) (W=200 -> xi ~ 23)
    while committed cells (pins and their persistently-driven shells)
    keep their full rate.
  * GEOMETRY FROZEN in v1 (labels dynamical, Pachner moves off): isolates
    the field physics; conservative field transport across moves is the
    documented v2 extension.

USAGE:
  python v6_capacity_run.py --resume clo_b1.0_20k_frozen.json \
         --pins-per-level 25 --sweeps 4000 --out cap_run
  python v6_capacity_run.py --analyze cap_run
"""
from __future__ import annotations
import argparse
import csv
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

from v6_run_lib import load_checkpoint
from v6_closure_run import (ETA_STAR, build_energy_table, tet_of, tris_of,
                            FaceLabels, _heatbath_pass)

G_SHARE_EFF = 7.4198          # vacuum budget anchor (paper App. B; verified)

# commitment ladder: frozen label configs with 0, 1, 2, 3, 6 colliding pairs.
# Level 0 (frozen but perfectly closed) is the DRESSING BASELINE: any frozen
# boundary slightly elevates neighbor failure regardless of its own
# commitment, contributing a level-independent flux offset. Q(m) - Q(0) is
# the commitment-scaling maintenance flux the paper's postulate concerns.
LADDER = {
    0: (3, 4, 5, 6),          # injective, K^2 minimal -- placebo rung
    1: (3, 3, 4, 5),          # one equal pair
    2: (3, 3, 4, 4),          # two equal pairs
    3: (3, 3, 3, 4),          # three equal pairs (triple)
    6: (3, 3, 3, 3),          # all equal -- the stage-3 'fail' defect
}


def build_slices(T):
    """Per-slice cell complexes: cells, face-adjacency, BFS-ready."""
    by_slice = defaultdict(set)
    for vs in T.pent.values():
        tet = tet_of(T, vs)
        if tet is not None:
            by_slice[T.vtime[next(iter(tet))]].add(tet)
    slices = {}
    for t, tets in by_slice.items():
        tets = sorted(tets, key=sorted)
        idx = {tet: i for i, tet in enumerate(tets)}
        tmap = defaultdict(list)
        for tet in tets:
            for tri in tris_of(tet):
                tmap[tri].append(idx[tet])
        adj = [[] for _ in tets]
        for cells in tmap.values():
            if len(cells) == 2:
                a, b = cells
                adj[a].append(b)
                adj[b].append(a)
        slices[t] = {"tets": tets, "idx": idx, "adj": adj,
                     "tet_pids": {tet: {0} for tet in tets}}
    return slices


def n_coll_vec(sl, labels, Etab):
    out = np.zeros(len(sl["tets"]))
    for i, tet in enumerate(sl["tets"]):
        l = [labels.lab[tri] for tri in tris_of(tet)]
        out[i] = sum(1 for a in range(4) for b in range(a + 1, 4)
                     if l[a] == l[b])
    return out


def bfs_dist(adj, src, cap=10**9):
    n = len(adj)
    dist = np.full(n, -1)
    dist[src] = 0
    frontier = [src]
    d = 0
    while frontier and d < cap:
        d += 1
        nxt = []
        for u in frontier:
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = d
                    nxt.append(v)
        frontier = nxt
    return dist


CSV_FIELDS = ["sweep", "level", "pin_id", "shell", "n_cells", "mean_f",
              "Q_pin", "far_f"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analyze", default=None, metavar="PREFIX",
                    help="analyze a finished run's CSV and exit")
    ap.add_argument("--resume", default=None,
                    help="thermalized closure checkpoint (geometry + labels)")
    ap.add_argument("--pins-per-level", type=int, default=25,
                    help="defects per commitment level (levels: 1,2,3,6)")
    ap.add_argument("--pin-sep", type=int, default=6)
    ap.add_argument("--rmax", type=int, default=5)
    ap.add_argument("--beta-closure", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=ETA_STAR)
    ap.add_argument("--lambda-inj", type=float, default=3.0)
    ap.add_argument("--D", type=float, default=0.2,
                    help="transport coefficient (stability: <= 0.25)")
    ap.add_argument("--kappa", type=float, default=0.05,
                    help="absorption per colliding pair per sweep")
    ap.add_argument("--absorb", choices=["excess", "total"], default="excess",
                    help="v1.1 default 'excess': absorb only failure ABOVE the "
                         "vacuum mean (measured during therm), so the vacuum "
                         "is absorption-free and the field massless -- required "
                         "for the T3/T4 Coulomb gates and for a proportional "
                         "M1. 'total' reproduces v1: vacuum turnover absorbs "
                         "too, giving intrinsic Yukawa screening at "
                         "xi = sqrt(D/(kappa*<n_coll>_vac)) ~ 5-7 steps (this "
                         "quantitatively explains the v1 run's observed 4-5 "
                         "step range and its M1 offset b~0.4).")
    ap.add_argument("--persist", type=int, default=200,
                    help="persistence window (EWMA sweeps) for the failure "
                         "average that drives absorption in excess mode; "
                         "1 = instantaneous (v1.1 behavior, intrinsic "
                         "xi ~ 7.5 -- see the module docstring)")
    ap.add_argument("--commit-margin", type=float, default=0.0,
                    help="optional threshold above nbar before a cell's "
                         "persistent failure draws maintenance (strict "
                         "commitment reading; 0 = soft persistence only)")
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--therm", type=int, default=400,
                    help="sweeps before measurements (field equilibration)")
    ap.add_argument("--post-therm", type=int, default=None,
                    help="field burn-in sweeps between the end of therm and "
                         "the first measurement (default: therm in --absorb "
                         "excess mode, 0 in total mode). In excess mode "
                         "absorption is OFF during therm -- the field stays "
                         "exactly uniform at the vacuum anchor while the "
                         "collision baseline accumulates -- and the massless "
                         "profile then needs ~rmax^2/D sweeps to form before "
                         "it is measured; measuring the buildup transient "
                         "would bias the T3/T4 fits toward screening.")
    ap.add_argument("--measure-every", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="cap_run")
    args = ap.parse_args()

    if args.analyze:
        analyze(args.analyze)
        return
    if not args.resume:
        sys.exit("--resume required (thermalized closure checkpoint)")

    Etab = build_energy_table(args.eta, args.lambda_inj)
    rng = np.random.default_rng(args.seed)
    T, _, extra = load_checkpoint(args.resume)
    labels = FaceLabels(rng)
    slices = build_slices(T)
    # ensure every triangle has a label; restore thermalized ones if present
    for sl in slices.values():
        for tet in sl["tets"]:
            labels.ensure(tris_of(tet))
    if isinstance(extra, dict) and extra.get("kind") == "closure":
        labels.load(extra)
        for sl in slices.values():
            for tet in sl["tets"]:
                labels.ensure(tris_of(tet))

    # --- pin the commitment ladder into the biggest slices ------------------
    levels = sorted(LADDER)
    order = sorted(slices, key=lambda t: -len(slices[t]["tets"]))
    pins = []                                # (level, t, cell_index)
    li = 0
    for t in order:
        sl = slices[t]
        adj = sl["adj"]
        taken = []
        cand = rng.permutation(len(sl["tets"]))
        for c in cand:
            if len(taken) >= max(1, (args.pins_per_level * len(levels))
                                 // max(1, len(order))) + 1:
                break
            if all(bfs_dist(adj, c, args.pin_sep)[o] < 0
                   or bfs_dist(adj, c, args.pin_sep)[o] >= args.pin_sep
                   for o in taken):
                taken.append(int(c))
        for c in taken:
            lvl = levels[li % len(levels)]
            li += 1
            pins.append((lvl, t, c))
            tet = sl["tets"][c]
            for tri, l in zip(tris_of(tet), LADDER[lvl]):
                labels.lab[tri] = int(l)
                labels.pinned.add(tri)
    per_lvl = {lvl: sum(1 for l, _, _ in pins if l == lvl) for lvl in levels}
    print(f"# [stage4] pinned ladder {per_lvl} across {len(order)} slices "
          f"(sep >= {args.pin_sep})", flush=True)

    post = args.post_therm if args.post_therm is not None else \
        (args.therm if args.absorb == "excess" else 0)
    meas_start = args.therm + post
    print(f"# [stage4] measurement window: sweeps {meas_start + 1}.."
          f"{args.sweeps}  (therm {args.therm} + field burn-in {post}; "
          f"persist {args.persist}, margin {args.commit_margin})",
          flush=True)
    if args.absorb == "excess" and meas_start < 3 * args.persist:
        print(f"# !! therm+post-therm ({meas_start}) < 3*persist "
              f"({3 * args.persist}): the persistence average may not have "
              f"converged when measurements start -- raise --therm/"
              f"--post-therm or lower --persist", flush=True)

    # --- field init (vacuum anchor T2) --------------------------------------
    f = {t: np.full(len(sl["tets"]), G_SHARE_EFF) for t, sl in slices.items()}
    tot0 = sum(v.sum() for v in f.values())
    # persistent failure average (EWMA), builds from sweep 1
    A = {t: np.zeros(len(sl["tets"])) for t, sl in slices.items()}

    # precompute BFS shells around each pin
    shells = []
    for lvl, t, c in pins:
        d = bfs_dist(slices[t]["adj"], c)
        shells.append((lvl, t, c, d))

    csv_f = open(args.out + ".csv", "w", newline="")
    writer = csv.DictWriter(csv_f, fieldnames=CSV_FIELDS)
    writer.writeheader()
    Q_acc = defaultdict(float)               # pin_id -> absorbed near pin
    Q_n = 0

    # vacuum baseline failure rate (for --absorb excess): accumulated over
    # the second half of therm, over cells far from every pin
    pin_far = {}
    for lvl, t, c, d in shells:
        pin_far.setdefault(t, np.ones(len(slices[t]["tets"]), bool))
        pin_far[t] &= ~((d >= 0) & (d <= args.rmax))
    nbar_acc, nbar_n = 0.0, 0

    for sw in range(1, args.sweeps + 1):
        # (i) label heat-bath per slice (beta-consistent; pins frozen)
        for t, sl in slices.items():
            _heatbath_pass(labels, sl["tet_pids"], Etab,
                           args.beta_closure, rng)
        # vacuum-baseline calibration during the second half of therm
        if args.therm // 2 < sw <= args.therm:
            for t, sl in slices.items():
                far = pin_far.get(t)
                ncol = n_coll_vec(sl, labels, Etab)
                sel = ncol[far] if far is not None else ncol
                nbar_acc += float(sel.sum()); nbar_n += len(sel)
            if sw == args.therm:
                nbar = nbar_acc / max(1, nbar_n)
                print(f"# [stage4] vacuum baseline <n_coll> = {nbar:.4f}  "
                      f"(absorb mode: {args.absorb})", flush=True)
        nbar_now = (nbar_acc / max(1, nbar_n)) if nbar_n else 0.0
        # (ii) transport + (iii) failure-proportional absorption, recycled
        for t, sl in slices.items():
            adj = sl["adj"]
            ff = f[t]
            flux = np.zeros_like(ff)
            for i, nbrs in enumerate(adj):
                for j in nbrs:
                    flux[i] += ff[j] - ff[i]
            ff = ff + args.D * flux
            ncol = n_coll_vec(sl, labels, Etab)
            A[t] += (ncol - A[t]) / args.persist
            if args.absorb == "excess":
                # no absorption during therm: the field stays exactly at the
                # uniform vacuum anchor while nbar accumulates (nbar is a
                # label-side statistic, independent of f), so the massless
                # profile builds from a clean initial condition. Rate is
                # driven by the PERSISTENT average A, not instantaneous
                # n_coll -- see the docstring (v1.2 rectifier fix).
                rate = (np.maximum(0.0, A[t] - nbar_now - args.commit_margin)
                        if sw > args.therm else np.zeros_like(ncol))
            else:
                rate = ncol
            absorb = args.kappa * rate * ff
            ff = ff - absorb + absorb.sum() / len(ff)
            f[t] = ff
        # T1 gate: exact conservation
        tot = sum(v.sum() for v in f.values())
        if abs(tot - tot0) > 1e-6 * tot0:
            raise RuntimeError(f"conservation violated at sweep {sw}: "
                               f"{tot0} -> {tot}")
        # measurements
        if sw > meas_start and sw % args.measure_every == 0:
            Q_n += 1
            for pid, (lvl, t, c, d) in enumerate(shells):
                sl = slices[t]
                # flux drawn by the pin + its immediate shell (emergent Q),
                # measured with the SAME rate law as the dynamics
                near = (d >= 0) & (d <= 1)
                if args.absorb == "excess":
                    rate_n = np.maximum(0.0, A[t][near] - nbar_now
                                        - args.commit_margin)
                else:
                    rate_n = n_coll_vec(sl, labels, Etab)[near]
                Q = float((args.kappa * rate_n * f[t][near]).sum())
                Q_acc[pid] += Q
                # same-slice far-field reference: removes the closed-slice
                # zero mode (recycled capacity elevates the far field) from
                # the deficit definition in the analyzer
                farm = pin_far.get(t)
                far_f = (round(float(f[t][farm].mean()), 6)
                         if farm is not None and farm.any() else "")
                for sh in range(1, args.rmax + 1):
                    m = d == sh
                    if m.any():
                        writer.writerow({
                            "sweep": sw, "level": lvl, "pin_id": pid,
                            "shell": sh, "n_cells": int(m.sum()),
                            "mean_f": round(float(f[t][m].mean()), 6),
                            "Q_pin": round(Q, 6), "far_f": far_f})
            csv_f.flush()
        if sw % 500 == 0:
            print(f"# sweep {sw}/{args.sweeps}  <f>="
                  f"{tot / sum(len(s['tets']) for s in slices.values()):.4f}"
                  f"  (conservation ok)", flush=True)

    csv_f.close()
    print(f"\n# STAGE 4 RUN DONE -> {args.out}.csv")
    print(f"# analyze:  python v6_capacity_run.py --analyze {args.out}")


def analyze(prefix):
    rows = list(csv.DictReader(open(prefix + ".csv")))
    if not rows:
        sys.exit("no rows")
    # deficit profile per level + emergent Q per level. Deficits are taken
    # against the SAME-SLICE far field when the CSV carries it (v1.2+;
    # removes the closed-slice zero mode), else against the vacuum anchor.
    def ref_of(r):
        v = r.get("far_f")
        return float(v) if v not in (None, "") else G_SHARE_EFF
    prof = defaultdict(list)                 # (level, shell) -> deficit
    Qs = defaultdict(list)                   # level -> Q per measurement
    for r in rows:
        prof[(int(r["level"]), int(r["shell"]))].append(
            ref_of(r) - float(r["mean_f"]))
        if int(r["shell"]) == 1:
            Qs[int(r["level"])].append(float(r["Q_pin"]))
    levels = sorted({int(r["level"]) for r in rows})
    shells = sorted({int(r["shell"]) for r in rows})
    print("shell DEFICITS (reference: same-slice far field when available; "
          "negative = elevation above it)")
    print(f"{'level':>6} " + "".join(f"{'d=%d' % s:>12}" for s in shells)
          + f" {'Q (emergent)':>14}")
    far = {}
    for lvl in levels:
        cells = []
        for sh in shells:
            v = np.array(prof[(lvl, sh)])
            cells.append(f"{v.mean():>12.4f}")
        q = np.array(Qs[lvl])
        far[lvl] = (np.array(prof[(lvl, shells[0])]).mean(),
                    q.mean())
        print(f"{lvl:>6} " + "".join(cells)
              + f" {q.mean():>10.4f} ±{q.std()/max(1,len(q))**0.5:.4f}")
    # M1: sink linearity  Q(level) ~ a*level + b, plus the proportional fit
    # (the paper's maintenance postulate wants Q proportional to commitment;
    # in --absorb excess mode b should be consistent with 0, in total mode
    # b is the vacuum-baseline absorption of the pin zone)
    ls = np.array(levels, float)
    qs = np.array([far[l][1] for l in levels])
    a, b = np.polyfit(ls, qs, 1)
    resid = float(np.abs(qs - (a * ls + b)).max() / max(qs.max(), 1e-12))
    a0 = float((qs * ls).sum() / (ls * ls).sum())
    resid0 = float(np.abs(qs - a0 * ls).max() / max(qs.max(), 1e-12))
    print(f"\nM1 sink linearity (affine):       Q = {a:.4f}*m + {b:.4f}  "
          f"(max rel. residual {resid:.3f})")
    print(f"M1 sink linearity (proportional): Q = {a0:.4f}*m           "
          f"(max rel. residual {resid0:.3f})")
    if 0 in far:
        b0 = far[0][1]
        m_pos = ls > 0
        qc = qs[m_pos] - b0
        lp = ls[m_pos]
        ac = float((qc * lp).sum() / (lp * lp).sum())
        rc = float(np.abs(qc - ac * lp).max() / max(abs(qc).max(), 1e-12))
        print(f"M1 DRESSING-SUBTRACTED (the claim): Q - Q(0) = {ac:.4f}*m  "
              f"(dressing baseline Q(0) = {b0:.4f}; max rel. residual {rc:.3f})")
    # M2 junction: shell-1 deficit per unit emergent Q (lattice-G analog).
    # Level 0 is the dressing baseline, not a mass -- excluded; both deficit
    # and flux are dressing-subtracted when the level-0 rung is present.
    b0 = far[0][1] if 0 in far else 0.0
    d0 = far[0][0] if 0 in far else 0.0
    lv = [l for l in levels if l > 0]
    d1 = np.array([far[l][0] - d0 for l in lv])
    qq = np.array([far[l][1] - b0 for l in lv])
    good = qq > 1e-9
    if good.any():
        Gl = (d1[good] / qq[good])
        print(f"M2 junction (dressing-subtracted shell-1 deficit per unit Q):"
              f" {Gl.mean():.4f} ± {Gl.std():.4f}  across levels "
              f"(constancy across the ladder = the lattice-G statement)")
    # T3/T4: profile-form and screening fits on the dressing-subtracted
    # deficit SHAPE (each level normalized to its own shell-1 deficit, then
    # pooled over levels > 0; second half of the measurements only)
    sw_all = sorted({int(r["sweep"]) for r in rows})
    lo = sw_all[len(sw_all) // 2]
    late = defaultdict(list)
    for r in rows:
        if int(r["sweep"]) >= lo:
            late[(int(r["level"]), int(r["shell"]))].append(
                ref_of(r) - float(r["mean_f"]))
    base = ({sh: float(np.mean(late[(0, sh)]))
             for sh in shells if late.get((0, sh))} if 0 in levels else {})
    if base:
        print("level-0 dressing profile (deficit; negative = elevation): "
              + "  ".join(f"d={sh}: {v:+.4f}" for sh, v in sorted(base.items())))
    shape = {}
    for sh in shells:
        vals = []
        for l in lv:
            v, v1 = late.get((l, sh)), late.get((l, 1))
            if v and v1:
                dsh = float(np.mean(v)) - base.get(sh, 0.0)
                d1 = float(np.mean(v1)) - base.get(1, 0.0)
                if d1 > 1e-12:
                    vals.append(dsh / d1)
        if vals:
            shape[sh] = (float(np.mean(vals)),
                         float(np.std(vals) / max(1, len(vals)) ** 0.5))
    print("\nT3/T4 deficit shape (dressing-subtracted, normalized to shell 1,"
          " pooled over levels > 0, late half):")
    for sh in sorted(shape):
        m, e = shape[sh]
        print(f"  d={sh}:  {m:+.4f} ± {e:.4f}")
    pos = [(sh, m) for sh, (m, _) in sorted(shape.items()) if m > 0]
    if len(pos) >= 3:
        x = np.array([p[0] for p in pos], float)
        y = np.array([p[1] for p in pos], float)
        co, cov = np.polyfit(np.log(x), np.log(y), 1, cov=True)
        p1, perr = float(co[0]), float(cov[0][0]) ** 0.5
        co2, cov2 = np.polyfit(x, np.log(y * x), 1, cov=True)
        s1, serr = float(co2[0]), float(cov2[0][0]) ** 0.5
        print(f"T3 profile form: deficit ~ d^({p1:.2f} ± {perr:.2f})   "
              f"(massless graph Green function on a 3-slice: ~ -1, with "
              f"flattening at the slice radius)")
        if abs(s1) < 2 * serr:
            bound = abs(s1) + 2 * serr
            print(f"T4 screening:    d ln(deficit*d)/dd = {s1:+.4f} ± {serr:.4f}"
                  f" -- consistent with MASSLESS (xi > {1 / bound:.1f} steps):"
                  f" gate PASS")
        elif s1 < 0:
            print(f"T4 screening:    d ln(deficit*d)/dd = {s1:+.4f} ± {serr:.4f}"
                  f" -- xi = {-1 / s1:.1f} steps: SCREENED (Route-A leak, or "
                  f"buildup transient if --post-therm was small): gate FAIL")
        else:
            print(f"T4 screening:    d ln(deficit*d)/dd = {s1:+.4f} ± {serr:.4f}"
                  f" -- positive (closed-slice zero mode / far-field surplus "
                  f"dominating): inspect the shape table")
        print("  caveats: the closed-slice zero mode elevates the far field, "
              "so negative-shape shells are excluded from the fits; quoted "
              "errors ignore sweep autocorrelation (thin with "
              "--measure-every if in doubt); usable range is capped by the "
              "slice radius -- check n_cells per shell in the CSV.")
    else:
        print("T3/T4: fewer than 3 shells with positive deficit -- range too "
              "short to fit (larger slices / larger N41 needed).")


if __name__ == "__main__":
    main()
