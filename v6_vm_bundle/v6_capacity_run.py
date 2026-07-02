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
    (iii) absorption dA_x = kappa * n_coll(x) * f_x, recycled uniformly
    within the slice the same sweep (steady state possible on a closed
    manifold).
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

# commitment ladder: frozen label configs with 1, 2, 3, 6 colliding pairs
LADDER = {
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
              "Q_pin"]


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
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--therm", type=int, default=400,
                    help="sweeps before measurements (field equilibration)")
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

    # --- field init (vacuum anchor T2) --------------------------------------
    f = {t: np.full(len(sl["tets"]), G_SHARE_EFF) for t, sl in slices.items()}
    tot0 = sum(v.sum() for v in f.values())

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

    for sw in range(1, args.sweeps + 1):
        # (i) label heat-bath per slice (beta-consistent; pins frozen)
        for t, sl in slices.items():
            _heatbath_pass(labels, sl["tet_pids"], Etab,
                           args.beta_closure, rng)
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
            absorb = args.kappa * ncol * ff
            ff = ff - absorb + absorb.sum() / len(ff)
            f[t] = ff
        # T1 gate: exact conservation
        tot = sum(v.sum() for v in f.values())
        if abs(tot - tot0) > 1e-6 * tot0:
            raise RuntimeError(f"conservation violated at sweep {sw}: "
                               f"{tot0} -> {tot}")
        # measurements
        if sw > args.therm and sw % args.measure_every == 0:
            Q_n += 1
            for pid, (lvl, t, c, d) in enumerate(shells):
                sl = slices[t]
                ncol = n_coll_vec(sl, labels, Etab)
                # flux drawn by the pin + its immediate shell (emergent Q)
                near = (d >= 0) & (d <= 1)
                Q = float((args.kappa * ncol[near] * f[t][near]).sum())
                Q_acc[pid] += Q
                for sh in range(1, args.rmax + 1):
                    m = d == sh
                    if m.any():
                        writer.writerow({
                            "sweep": sw, "level": lvl, "pin_id": pid,
                            "shell": sh, "n_cells": int(m.sum()),
                            "mean_f": round(float(f[t][m].mean()), 6),
                            "Q_pin": round(Q, 6)})
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
    # deficit profile per level + emergent Q per level
    prof = defaultdict(list)                 # (level, shell) -> mean_f
    Qs = defaultdict(list)                   # level -> Q per measurement
    for r in rows:
        prof[(int(r["level"]), int(r["shell"]))].append(float(r["mean_f"]))
        if int(r["shell"]) == 1:
            Qs[int(r["level"])].append(float(r["Q_pin"]))
    levels = sorted({int(r["level"]) for r in rows})
    shells = sorted({int(r["shell"]) for r in rows})
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
    # M1: sink linearity  Q(level) ~ a*level + b
    ls = np.array(levels, float)
    qs = np.array([far[l][1] for l in levels])
    a, b = np.polyfit(ls, qs, 1)
    resid = float(np.abs(qs - (a * ls + b)).max() / max(qs.max(), 1e-12))
    print(f"\nM1 sink linearity: Q = {a:.4f}*commitment + {b:.4f}   "
          f"(max rel. residual {resid:.3f})")
    # M2 junction: shell-1 deficit per unit emergent Q (lattice-G analog)
    d1 = np.array([G_SHARE_EFF - far[l][0] for l in levels])
    good = qs > 1e-9
    if good.any():
        Gl = (d1[good] / qs[good])
        print(f"M2 junction (shell-1 deficit per unit Q): "
              f"{Gl.mean():.4f} ± {Gl.std():.4f}  across levels "
              f"(constancy across the ladder = the lattice-G statement)")
    print("\nT3/T4 (profile form + screening fit) need the larger-volume "
          "run: fit ln(deficit) vs ln(d) / d across shells at 20k.")


if __name__ == "__main__":
    main()
