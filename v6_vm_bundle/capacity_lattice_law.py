#!/usr/bin/env python3
"""FROM THEORY TOWARD REGGE, step 4: the lattice capacity-conservation law
-- what the paper's finite-capacity budget implies on the lattice, and the
Coulomb/Yukawa dichotomy that decides whether the theory can have long-range
gravity at all. (Derivation in CAPACITY_CONSERVATION.md; this script checks
the anchors and demonstrates the dichotomy numerically on a REAL slice
geometry from a run checkpoint.)

WHAT THE PAPER FIXES (inputs to the derivation):
  P1  finite capacity budget per cell; the vacuum cell's sharing entropy is
      the admissibility-closed value g_share,eff = 7.4198 (verified below --
      this is the budget anchor, same number as the closure-model gate);
  P2  commitment: defects lock capacity; committed capacity is conserved
      except at commitment/release events (mass conservation);
  P3  defects need CONTINUOUS MAINTENANCE, linear in commitment (the paper's
      own dressing statement: maintenance is large and linear in the
      commitment while rest energy is the rare-return gap);
  P4  the refresh supplies each cell a fixed redraw bandwidth 1/tau*.

THE DICHOTOMY (the theorem-shaped result):
  Route A -- BUDGET ONLY. If free capacity is merely constrained per cell
  (sum of link shares = budget - commitment) and re-equilibrates by MaxEnt,
  the linear response to a point commitment solves (2z - L) dmu = -m/chi
  (z = 4 faces, L = slice Laplacian): an operator with NO small-k pole.
  Response is SCREENED with sub-cell range. Gravity would be exponentially
  cut off at the substrate scale -- no Newton, ever.

  Route B -- CONSERVED TRANSPORT + STEADY SINK. If total free capacity obeys
  a continuity law (changed only by face-fluxes and commit/release) and a
  defect draws a steady maintenance flux Q proportional to its commitment
  (P3), the steady state solves L drho = (Q/D)(delta_source - 1/N): the
  MASSLESS graph-Poisson equation. On a 3D slice its Green function falls
  like 1/r -- Newtonian form, with Phi ~ Q/r ~ (commitment)/r, i.e. the
  M-linearity of Newton's law inherited directly from P3.

  CONCLUSION: the paper's long-range gravity REQUIRES the conservation
  reading (Route B). Budget bookkeeping alone is not enough; the lattice
  law any propagation model must implement is the continuity equation

      f_x(t+dt) - f_x(t) = - sum_faces J_xy(t)  - commit_x + release_x ,
      J antisymmetric,  sum_x [f_x + m_x] = const (exact, machine precision)

  with vacuum anchor <f> = g_share,eff = 7.4198 per cell and sink strength
  linear in m_x. These are the fidelity gates (T1-T4 in the .md) for any
  future stage-4 'propagation' model, playing the role 7.4198 played for
  the closure model.

Usage:  python capacity_lattice_law.py [checkpoint.json]
        (defaults to a fresh S1xS3 slice if no checkpoint is given)
"""
from __future__ import annotations
import sys
from collections import defaultdict, Counter
from itertools import combinations, permutations

import numpy as np

from v6_closure_run import ETA_STAR, tet_of, tris_of


def anchor_entropy(eta=ETA_STAR):
    """The budget anchor: entropy of the admissibility ensemble (paper's
    g_share,eff = 7.4198)."""
    K2 = []
    for m in permutations(range(-3, 4), 4):
        S = sum(m); Sig2 = sum(x * x for x in m)
        K2.append(48.0 - (S * S - Sig2) / 3.0)
    K2 = np.array(K2)
    w = np.exp(-eta * K2)
    Z = w.sum()
    p = w / Z
    # over 1680 states (parity doubling halves probabilities):
    return float(-(p * np.log(p / 2.0)).sum())


def biggest_slice_graph(path):
    """Cell-adjacency graph of the largest time slice in a checkpoint."""
    from v6_run_lib import load_checkpoint
    T, _, _ = load_checkpoint(path)
    by_slice = defaultdict(set)
    for vs in T.pent.values():
        tet = tet_of(T, vs)
        if tet is not None:
            t = T.vtime[next(iter(tet))]
            by_slice[t].add(tet)
    t_big = max(by_slice, key=lambda t: len(by_slice[t]))
    tets = sorted(by_slice[t_big], key=sorted)
    idx = {tet: i for i, tet in enumerate(tets)}
    tmap = defaultdict(list)
    for tet in tets:
        for tri in tris_of(tet):
            tmap[tri].append(idx[tet])
    n = len(tets)
    A = np.zeros((n, n))
    for tri, cells in tmap.items():
        for a in cells:
            for b in cells:
                if a != b:
                    A[a, b] = 1.0
    return A, t_big


def dichotomy_demo(A):
    """Both response profiles on the same real slice graph."""
    n = A.shape[0]
    deg = A.sum(1)
    L = np.diag(deg) - A
    z = float(deg.mean())
    src = 0

    # BFS distances from the source
    dist = np.full(n, -1); dist[src] = 0
    frontier = [src]
    d = 0
    while frontier:
        d += 1
        nxt = []
        for u in frontier:
            for v in np.nonzero(A[u])[0]:
                if dist[v] < 0:
                    dist[v] = d; nxt.append(v)
        frontier = nxt

    # Route A: (2z - L) x = e_src   (budget-only MaxEnt response)
    xa = np.linalg.solve(2 * z * np.eye(n) - L, np.eye(n)[src])
    # Route B: L x = e_src - 1/n    (conserved transport, steady sink,
    #                                uniform recycling; pseudo-inverse)
    rhs = np.eye(n)[src] - 1.0 / n
    xb = np.linalg.lstsq(L, rhs, rcond=None)[0]
    xb -= xb.mean()

    print(f"\n[dichotomy on the real slice graph: n={n} cells, mean degree "
          f"{z:.1f}]")
    print(f"{'BFS dist':>8} {'cells':>6} {'Route A (budget-only)':>22} "
          f"{'Route B (conserved)':>21}")
    a0 = abs(xa[src]); b0 = abs(xb[src])
    for dd in range(0, int(dist.max()) + 1):
        m = dist == dd
        if not m.any():
            continue
        print(f"{dd:>8} {int(m.sum()):>6} {abs(xa[m]).mean() / a0:>22.2e} "
              f"{abs(xb[m]).mean() / b0:>21.2e}")
    print("\n  Route A collapses by orders of magnitude per step (screened,"
          " sub-cell range):\n  no long-range gravity is possible from budget"
          " bookkeeping alone.")
    print("  Route B decays like the massless graph-Coulomb Green function"
          " (~1/r on a 3D\n  slice, flattened here by the small closed"
          " manifold): the conservation law is\n  what carries strain to"
          " long range. This is the lattice law a propagation model\n  must"
          " implement.")


if __name__ == "__main__":
    g = anchor_entropy()
    print(f"[anchor] admissibility-ensemble entropy = {g:.4f} "
          f"(paper g_share,eff = 7.4198) -> per-cell vacuum budget")
    if len(sys.argv) > 1:
        A, t = biggest_slice_graph(sys.argv[1])
        print(f"[slice] using largest slice (t={t}) of {sys.argv[1]}")
    else:
        from v6_cdt import build_s1xs3
        import tempfile, json, os
        T = build_s1xs3(K=8)
        from v6_run_lib import save_checkpoint
        tmp = tempfile.mktemp(suffix=".json")
        save_checkpoint(T, tmp)
        A, t = biggest_slice_graph(tmp)
        os.unlink(tmp)
        print(f"[slice] no checkpoint given -- fresh S1xS3 slice (t={t})")
    dichotomy_demo(A)
