#!/usr/bin/env python3
"""STAGE 3 PREVIEW (free -- runs on existing stage-2 checkpoints): do
naturally-occurring transient closure failures sit in geometrically different
environments than well-closed cells?

This is the fluctuation-response version of the pinned-defect experiment
(v6_defect_run.py): instead of inserting a persistent source, use the ~8% of
cells that happen to be failing at the checkpointed instant and compare their
local environment against the closing cells' -- same physics by
fluctuation-dissipation logic, much weaker statistics (ONE label snapshot per
checkpoint, since checkpoints overwrite).

Per cell: failing = any label collision among its four faces.
Environment observables:
  mean_q    mean coordination of the cell's 6 slice edges (local curvature)
  ball2     number of tets within BFS distance 2 (local volume)
  nbrE      mean closure energy of the 4 face-neighbors (capacity strain)
Report: failing vs closing means with bootstrap errors.

Usage:  python failure_geometry_preview.py clo_b1.0_20k.json [more ckpts...]
"""
from __future__ import annotations
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

from v6_run_lib import load_checkpoint
from v6_closure_run import (ETA_STAR, build_energy_table, tet_of, tris_of)


def preview(path, eta=ETA_STAR, lam=3.0):
    T, meta, extra = load_checkpoint(path)
    if not (isinstance(extra, dict) and extra.get("kind") == "closure"):
        print(f"{path}: not a closure checkpoint (no labels) -- skipped")
        return
    lab = {frozenset((a, b, c)): l for a, b, c, l in extra["lab"]}
    E = build_energy_table(eta, lam)

    tets = set()
    for vs in T.pent.values():
        tet = tet_of(T, vs)
        if tet is not None:
            tets.add(tet)
    # adjacency + edge coordination over the slices
    tmap = defaultdict(list)
    ecount = defaultdict(int)
    for tet in tets:
        for tri in tris_of(tet):
            tmap[tri].append(tet)
        for e in combinations(sorted(tet), 2):
            ecount[e] += 1
    adj = defaultdict(set)
    for tri, ts in tmap.items():
        for a in ts:
            for b in ts:
                if a != b:
                    adj[a].add(b)

    def cell_E(tet):
        l = [lab[t] for t in tris_of(tet)]
        return float(E[l[0], l[1], l[2], l[3]]), len(set(l)) < 4

    groups = {True: defaultdict(list), False: defaultdict(list)}
    missing = 0
    for tet in tets:
        try:
            _, failing = cell_E(tet)
        except KeyError:
            missing += 1
            continue
        g = groups[failing]
        g["mean_q"].append(np.mean([ecount[e] for e in
                                    combinations(sorted(tet), 2)]))
        seen = {tet} | adj[tet]
        ball2 = set(seen)
        for u in adj[tet]:
            ball2 |= adj[u]
        g["ball2"].append(len(ball2) - 1)
        nbrEs = []
        for u in adj[tet]:
            try:
                ev, _ = cell_E(u)
                nbrEs.append(ev)
            except KeyError:
                pass
        if nbrEs:
            g["nbrE"].append(np.mean(nbrEs))

    nf, nc = len(groups[True]["mean_q"]), len(groups[False]["mean_q"])
    print(f"\n== {path} ==")
    print(f"  cells: {len(tets)}  failing: {nf} ({nf / max(1, nf + nc):.3f})"
          f"  (label-less skipped: {missing})")
    rng = np.random.default_rng(0)
    print(f"  {'observable':>10} {'failing':>18} {'closing':>18} {'diff':>18}")
    for key in ("mean_q", "ball2", "nbrE"):
        a = np.array(groups[True][key], float)
        b = np.array(groups[False][key], float)
        if len(a) < 5 or len(b) < 5:
            continue
        boots = [a[rng.integers(0, len(a), len(a))].mean()
                 - b[rng.integers(0, len(b), len(b))].mean()
                 for _ in range(400)]
        d, se = float(np.mean(boots)), float(np.std(boots))
        print(f"  {key:>10} {a.mean():>12.4f}       {b.mean():>12.4f}      "
              f"{d:>+9.4f} ± {se:.4f}"
              f"   {'<-- separated' if abs(d) > 2 * se else ''}")
    print("  (one snapshot -- preview-grade statistics; the pinned-defect "
          "run is the real measurement)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        preview(p)
