#!/usr/bin/env python3
"""FROM THEORY TOWARD REGGE, step 5: integrate out the capacity field --
what geometric action does the BACK-COUPLED conserved field induce, exactly?

THE QUESTION (raised by the observation that in this theory vacuum and
matter are one medium): the earlier no-gos closed the LABEL sector -- local
interactions, sub-lattice correlation lengths, structurally unable to order
geometry. The stage-4 conserved capacity field is a NEW degree of freedom
those proofs never covered: long-ranged by construction (massless because
conserved). If geometry feels the field (v2 back-coupling), what action does
the field induce on geometry?

THE ANSWER IS EXACT. A conserved field with gradient stiffness at
temperature T (fluctuation-dissipation: noisy transport; the v1 field is
deterministic and needs the noise term for this to apply) is Gaussian:

    Z_field(g) = integral df exp(-(kappa/2T) f L(g) f)
               = const^(N-1) * (det' L(g))^(-1/2)         (zero mode = the
                                                            conserved total)
    =>  Delta S_eff(g) = -log Z = (n_f/2) * log det' L(g)  + volume terms,

per field component n_f (the theory's capacity medium is 7-channel, so
n_f = 7 is the natural value). By Kirchhoff's theorem det' L = N * (number
of spanning trees), so equivalently: THE FIELD ENTROPICALLY FAVORS
GEOMETRIES WITH SMALL SPECTRAL DETERMINANTS -- rich in soft long-wavelength
modes -- i.e. LARGE, EXTENDED, SMOOTH slices, and it PENALIZES gapped
expander-like (crumpled / branched) geometry.

WHY THIS EVADES THE NO-GOS: log det' L is NOT a sum of local terms. Its
small-eigenvalue part is a functional of the slice's GLOBAL connectivity --
exactly the nonlocality the label sector provably lacked. This is the
standard mechanism by which matter fields shift dynamical-triangulation
phase diagrams, here arising from the theory's own medium.

WHAT THIS SCRIPT COMPUTES (all on 4-regular graphs, matched degree):
  phi(g) = (1/2N) log det' L(g)   -- induced action per cell -- for
  (a) diamond-cubic lattice (proxy for extended smooth 3-geometry: large
      diameter, soft modes),
  (b) random 4-regular graphs (proxy for crumpled/branched: expander,
      spectral gap),
  (c) REAL slice graphs from our run checkpoints (the engine's current,
      collapsed-phase slices),
and reports the per-cell differential Delta(phi) between crumpled and
extended -- the entropic force per cell per field component -- and the
channel count n_f at which the total reaches the k0 ~ 2.2 scale that
sustains phase C.

READING: if n_f ~ 7 puts the field's phase-tipping force at O(bare
couplings), the theory-native route to the de Sitter phase (geometry +
its own medium, v2 back-coupling) is quantitatively live, and 'the vacuum
works when the medium is present' becomes a testable prediction rather
than a hope. If the differential is tiny, the field is phase-irrelevant
and the no-go extends to the field sector too.

Usage:  python field_induced_action.py [checkpoint.json ...]
"""
from __future__ import annotations
import sys
from collections import defaultdict

import numpy as np


def phi_per_cell(A):
    """(1/2N) log det' L for adjacency matrix A (dense, small graphs)."""
    deg = A.sum(1)
    L = np.diag(deg) - A
    ev = np.linalg.eigvalsh(L)
    ev = ev[1:]                              # drop the conserved zero mode
    ev = np.clip(ev, 1e-12, None)
    return float(np.log(ev).sum() / (2 * len(A)))


def diamond_cubic(n):
    """Periodic diamond-cubic lattice (4-regular): 2 sites per fcc cell.
    n = cells per axis -> N = 2 n^3 sites."""
    N = 2 * n ** 3
    def idx(s, x, y, z):
        return ((x % n) * n + (y % n)) * n + (z % n) + s * n ** 3
    A = np.zeros((N, N))
    # site 0 at (0,0,0), site 1 at (1/4,1/4,1/4); bonds: 1-site connects to
    # 0-sites of its own cell and the three neighbor cells offset by -x,-y,-z
    for x in range(n):
        for y in range(n):
            for z in range(n):
                a = idx(1, x, y, z)
                for b in (idx(0, x, y, z), idx(0, x + 1, y, z),
                          idx(0, x, y + 1, z), idx(0, x, y, z + 1)):
                    A[a, b] = A[b, a] = 1
    return A


def random_regular(N, d=4, seed=0):
    """Random d-regular graph via pairing (retry until simple)."""
    rng = np.random.default_rng(seed)
    for _ in range(200):
        stubs = np.repeat(np.arange(N), d)
        rng.shuffle(stubs)
        A = np.zeros((N, N))
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or A[a, b]:
                ok = False
                break
            A[a, b] = A[b, a] = 1
        if ok:
            return A
    raise RuntimeError("no simple regular graph found")


def slice_graphs(path, max_slices=3):
    from v6_run_lib import load_checkpoint
    from v6_closure_run import tet_of, tris_of
    T, _, _ = load_checkpoint(path)
    by_slice = defaultdict(set)
    for vs in T.pent.values():
        tet = tet_of(T, vs)
        if tet is not None:
            by_slice[T.vtime[next(iter(tet))]].add(tet)
    out = []
    for t in sorted(by_slice, key=lambda t: -len(by_slice[t]))[:max_slices]:
        tets = sorted(by_slice[t], key=sorted)
        idx = {tet: i for i, tet in enumerate(tets)}
        tmap = defaultdict(list)
        for tet in tets:
            for tri in tris_of(tet):
                tmap[tri].append(idx[tet])
        A = np.zeros((len(tets), len(tets)))
        for cells in tmap.values():
            if len(cells) == 2:
                a, b = cells
                A[a, b] = A[b, a] = 1
        out.append((t, A))
    return out


if __name__ == "__main__":
    print(f"{'graph':<38} {'N':>6} {'phi = (1/2N) log det`L':>24}")
    rows = {}
    for n in (4, 5, 6):
        A = diamond_cubic(n)
        p = phi_per_cell(A)
        rows.setdefault("extended", []).append(p)
        print(f"{'diamond-cubic (extended smooth)':<38} {len(A):>6} {p:>24.4f}")
    for N, seed in ((128, 0), (250, 1), (432, 2)):
        A = random_regular(N, 4, seed)
        p = phi_per_cell(A)
        rows.setdefault("crumpled", []).append(p)
        print(f"{'random 4-regular (crumpled proxy)':<38} {N:>6} {p:>24.4f}")
    for path in sys.argv[1:]:
        for t, A in slice_graphs(path):
            if len(A) < 30:
                continue
            p = phi_per_cell(A)
            rows.setdefault("engine slices", []).append(p)
            print(f"{'REAL slice t=%d  %s' % (t, path.split('/')[-1]):<38} "
                  f"{len(A):>6} {p:>24.4f}")

    ext = np.mean(rows["extended"])
    cru = np.mean(rows["crumpled"])
    d = cru - ext
    print(f"\nper-cell entropic force (crumpled - extended) per field "
          f"component: {d:+.4f}")
    print(f"  (positive = the field PENALIZES crumpled geometry, as the "
          f"soft-mode argument predicts)")
    if d > 0:
        for target in (0.5, 2.2):
            print(f"  channels needed for total ~ {target}: "
                  f"n_f ~ {target / d:.1f}")
        print(f"  the theory's capacity medium is 7-channel: n_f = 7 gives "
              f"a phase-tipping force of ~ {7 * d:.2f} per cell")
    if "engine slices" in rows:
        eng = np.mean(rows["engine slices"])
        print(f"\n  engine's current slices sit at phi = {eng:.4f} "
              f"(extended = {ext:.4f}, crumpled = {cru:.4f}): the field "
              f"would push them by {(eng - ext):+.4f}/cell toward extended.")
