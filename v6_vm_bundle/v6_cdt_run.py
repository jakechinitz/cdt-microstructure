#!/usr/bin/env python3
"""v6 Metropolis driver -- phase 1: demonstrate the engine evolves and provably
stays a valid CDT 4-manifold.

Design notes (kept deliberately simple/readable):
  * Moves are applied in place; a REJECTED move is undone by applying its
    verified INVERSE (we never snapshot the whole state). Each move returns the
    handle its inverse needs (the new edge / new vertex / dual hinge).
  * Action = bare Regge + quadratic volume penalty on N_41 (recomputed each
    step -- O(N), fine at this scale; an O(1) delta + a bag proposer are the
    next optimisation).
  * Every `verify_every` sweeps we run the full independent verify(): the chain
    must remain a closed 4-manifold with chi=0 and S^3 vertex links. This is the
    robustness guarantee v5 could not give.

NOT yet here (next steps): exact detailed-balance Hastings factors, the
(4,6)/(6,4) move (full ergodicity), strictness tuning, the bag proposer, and the
phase-C (d_s->4) validation.
"""
from __future__ import annotations
import numpy as np
from itertools import combinations

from v6_cdt import build_s1xs3, verify, regge_action
from v6_cdt_moves import (
    apply_24, apply_42, apply_28, apply_82, apply_33, apply_46, apply_64,
    strictness_ok, causal_slice_ok,
)

MOVES = ["(2,4)", "(4,2)", "(2,8)", "(8,2)", "(3,3)", "(4,6)", "(6,4)"]


def _action(T, k0, D, k4, N_bar, eps):
    S = regge_action(T, k0, D, k4)
    if N_bar is not None and eps > 0:
        n41 = T.type_counts()[0]
        S += eps * (n41 - N_bar) ** 2
    return S


def hastings_log(mt, N4, N0):
    """Detailed-balance proposal ratio log[q(x'->x)/q(x->x')] for the proposal
    scheme in propose_and_apply (pick a uniform pentachoron + a uniform
    sub-simplex; (8,2) picks a uniform vertex). Derived from the move
    multiplicities; N4,N0 are the counts BEFORE the move. (3,3) is volume- and
    proposal-symmetric so its factor is 0."""
    from math import log
    if mt == "(2,4)":
        return log(N4) - log(N4 + 2)
    if mt == "(4,2)":
        return log(N4) - log(N4 - 2)
    if mt == "(2,8)":
        return log(5 * N4) - log(2 * (N0 + 1))
    if mt == "(8,2)":
        return log(2 * N0) - log(5 * (N4 - 6))
    if mt == "(4,6)":
        return log(3 * N4) - log(2 * (N4 + 2))
    if mt == "(6,4)":
        return log(2 * N4) - log(3 * (N4 - 2))
    return 0.0  # (3,3)


def propose_and_apply(T, rng, strictness=2, causal=True, protect=None):
    """Pick a uniform move type + a random target, apply it, and return
    (move_label, ok, undo). A move that would create a degenerate (non-strict)
    configuration is undone and reported as ok=False -- so the chain stays on the
    strict (non-branched) ensemble (detailed balance is unaffected: the proposal
    distribution is unchanged; strictness only forbids entering degenerate states).

    `causal=True` (default, 2026-07 fix) additionally rejects moves that break
    the CDT foliation (spatial slices must remain closed simplicial 3-manifolds
    with every spatial tetrahedron an interface between one future and one past
    pentachoron -- see causal_slice_ok). The pre-fix ensemble (causal=False)
    develops slice defects through (2,4)/(3,3) moves and is NOT the standard
    CDT ensemble of the literature.

    `protect` (optional set of vertex ids): veto (8,2) proposals on these
    vertices BEFORE applying. Needed by the stage-3 defect runs: an (8,2)
    that is applied and then undone (by ANY downstream rejection) re-creates
    the vertex with a fresh id, which silently unkeys vertex-identified
    structure like pinned defect carriers. Vetoing pre-apply is the only
    id-stable protection; it restricts the ensemble exactly like the other
    filters and must be identical across compared arms."""
    mt = MOVES[rng.integers(len(MOVES))]
    P = T.pent_bag.pick(rng)        # O(1) random pentachoron
    vs = T.pent[P]
    news = None
    undo = None

    if mt == "(2,4)":
        face = tuple(sorted(_rand_subset(vs, 4, rng)))
        ok, news = apply_24(T, P, face)
        if ok:
            e, f = sorted(set(T.pent[news[0]]) & set(T.pent[news[1]]) &
                          set(T.pent[news[2]]) & set(T.pent[news[3]]))
            undo = lambda: apply_42(T, e, f)
    elif mt == "(4,2)":
        e, f = sorted(_rand_subset(vs, 2, rng))
        ok, news = apply_42(T, e, f)
        if ok:
            Pn = news[0]
            abcd = tuple(sorted(set(T.pent[Pn]) - (set(T.pent[Pn]) - set(T.pent[news[1]]))))
            undo = lambda: apply_24(T, Pn, abcd)
    elif mt == "(2,8)":
        face = tuple(sorted(_rand_subset(vs, 4, rng)))
        ok = len({T.vtime[v] for v in face}) == 1
        if ok:
            ok, news, x = apply_28(T, P, face)
            if ok:
                undo = lambda: apply_82(T, x)
    elif mt == "(8,2)":
        x = T.vert_bag.pick(rng)        # O(1) random vertex
        if protect and x in protect:
            return mt, False, None      # protected carrier vertex (stage 3)
        ok, news = apply_82(T, x)
        if ok:
            Pn = news[0]
            abcd = tuple(sorted(set(T.pent[Pn]) - (set(T.pent[Pn]) - set(T.pent[news[1]]))))
            undo = lambda: apply_28(T, Pn, abcd)
    elif mt == "(3,3)":
        tri = tuple(sorted(_rand_subset(vs, 3, rng)))
        ok, news, dual = apply_33(T, P, tri)
        if ok:
            undo = lambda: apply_33(T, news[0], dual)
    elif mt == "(4,6)":
        tri = tuple(sorted(_rand_subset(vs, 3, rng)))
        ok = len({T.vtime[v] for v in tri}) == 1
        if ok:
            ok, news, ab = apply_46(T, P, tri)
            if ok:
                a, b = ab
                undo = lambda: apply_64(T, a, b)
    elif mt == "(6,4)":
        e = tuple(sorted(_rand_subset(vs, 2, rng)))
        ok = T.vtime[e[0]] == T.vtime[e[1]]
        if ok:
            ok, news, tri = apply_64(T, e[0], e[1])
            if ok:
                Pn = news[0]
                undo = lambda: apply_46(T, Pn, tuple(sorted(set(tri))))
    else:
        ok = False

    if not ok:
        return mt, False, None
    if not strictness_ok(T, news, strictness):
        undo()                       # would create a degenerate config -> reject
        return mt, False, None
    if causal and not causal_slice_ok(T, news):
        undo()                       # would break the foliation -> reject
        return mt, False, None
    return mt, True, undo


def _rand_subset(vs, k, rng):
    idx = rng.choice(len(vs), size=k, replace=False)
    return [vs[i] for i in idx]


def run(K=10, sweeps=60, k0=2.5, Delta=0.6, k4=1.0, N_bar=None, eps=1e-2,
        beta=1.0, seed=0, verify_every=10, verbose=True, strictness=2):
    rng = np.random.default_rng(seed)
    T = build_s1xs3(K=K)
    if N_bar is None:
        N_bar = T.type_counts()[0] * 3   # grow to ~3x the init N_41
    S = _action(T, k0, Delta, k4, N_bar, eps)
    acc = {m: 0 for m in MOVES}
    prop = {m: 0 for m in MOVES}
    ok0, _ = verify(T)
    if verbose:
        print(f"init: N4={T.n_pent()} N41={T.type_counts()[0]} verify={'PASS' if ok0 else 'FAIL'} "
              f"(target N41={N_bar})")

    for sw in range(1, sweeps + 1):
        for _ in range(max(1, T.n_pent())):
            N4_b, N0_b = T.n_pent(), len(T.vinc)
            mt, ok, undo = propose_and_apply(T, rng, strictness)
            prop[mt] += 1
            if not ok:
                continue
            S_new = _action(T, k0, Delta, k4, N_bar, eps)
            H = hastings_log(mt, N4_b, N0_b)
            if rng.random() < np.exp(min(50.0, -beta * (S_new - S) + H)):
                S = S_new
                acc[mt] += 1
            else:
                res = undo()                      # apply the inverse move
                u_ok = res[0] if isinstance(res, tuple) else res
                if not (u_ok and T.n_pent() == N4_b):
                    raise RuntimeError(
                        f"undo of {mt} FAILED (N4 {N4_b}->{T.n_pent()}) -- "
                        f"a move/inverse bug; refusing to continue with corrupt state")
        if sw % verify_every == 0:
            ok, rep = verify(T)
            if verbose:
                print(f"sweep {sw:>3}: N4={T.n_pent():>4} N41={T.type_counts()[0]:>4} "
                      f"verify={'PASS' if ok else 'FAIL'} chi={rep['euler_char']} "
                      f"links={rep['link_failures']} types_other={rep['type_counts'][2]}")
            if not ok:
                print("  !!! manifold validity BROKEN -- stopping")
                return T, False
    if verbose:
        print("accepts by move:", {m: f"{acc[m]}/{prop[m]}" for m in MOVES})
    return T, True


if __name__ == "__main__":
    print("=" * 72)
    print("  v6 Metropolis driver -- does it run and STAY a valid 4-manifold?")
    print("=" * 72)
    T, ok = run(K=10, sweeps=60, verify_every=10, seed=1)
    print(f"\nFINAL: valid-throughout={ok}  N4={T.n_pent()}")
