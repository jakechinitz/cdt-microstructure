#!/usr/bin/env python3
"""FROM THEORY TOWARD REGGE, step 3 (layer 2): the Many-Pasts inverse problem.

Steps 1-2 established a two-branch no-go: neither the static closure weight
(induced curvature coupling ~1% of need, link ratio r0 = |lam2/lam1| = 0.12)
nor Markovian refresh (geometry-blind NESS) supplies vacuum stiffness. The
one remaining ingredient is Postulate III's HISTORY WEIGHTING. This script
does not invent a Many-Pasts lattice dynamics; it computes what any such
weighting WOULD HAVE TO DO, in two parts:

PART A (model-independent requirement). The curvature-sector signal is the
q-dependence of the ring free energy F(q) (q = coordination = discrete
curvature). Whatever the microscopic origin, an effective label ensemble is
characterized by its link transfer matrix; amplifying the interaction
(effective inverse temperature beta_eff on the SAME closure energies) maps
out how much link correlation is needed. We report, vs beta_eff: the link
ratio r, the correlation length xi = -1/ln r, and the ring curvature signal
(per-edge linear coefficient c0 and the nonlinear remainder). REQUIREMENT =
the beta_eff / r / xi where the signal reaches O(1) (thresholds 0.5 and 2.2,
the bare k0 sustaining phase C).

PART B (a faithful-in-spirit realization, solved exactly). Many-Pasts
weights histories, not states -- mathematically a PATH-MEASURE TILT. The
closest closure-flavored reading: tilt the label dynamics by the
time-integrated closure quality, exp(-s * int V dt), V = sum of link
energies -- 'histories that maintain closure count more'. For a reversible
single-face heat-bath on a ring of q faces this is exactly solvable: the
tilted process's equal-time measure is psi_s^2 where psi_s is the leading
eigenvector of the symmetrized tilted generator H_s = H - s V (Doob
transform). We compute, vs tilting strength s:
    * the equal-time link correlation and its decay along the ring
      (does tilting actually amplify SPATIAL correlations? -> yes: the
      measure interpolates toward a ground-state-like object)
    * the COST: chi(s) = per-face rarity rate of the selected histories
      (how atypical the Many-Pasts ensemble must be, per cell per refresh
      time, to buy the correlation).

CAUSALITY NOTE (why this evades the light-cone objection): a causal local
dynamics with memory horizon tau_mem cannot sustain equal-time correlations
beyond ~tau_mem links. Path tilting conditions on histories globally --
including retrospectively -- which is precisely the acausal-in-form,
operationally-safe structure Many-Pasts claims. The requirement computed
here is therefore also a statement of HOW acausal-in-form the weighting
must be: correlation length xi_req demands history conditioning that no
xi_req-generation-deep causal kernel could mimic.

Usage:  python many_pasts_requirements.py [--q 6] [--fast]
"""
from __future__ import annotations
import argparse
from itertools import product
import numpy as np

from v6_closure_run import ETA_STAR, build_energy_table


def link_matrix(beta_eff, eta=ETA_STAR, lam=3.0):
    """Step-1 transfer matrix at effective inverse temperature beta_eff."""
    E = build_energy_table(eta, lam)
    T = np.zeros((7, 7))
    for m in range(7):
        for mp in range(7):
            T[m, mp] = np.exp(-beta_eff * E[m, mp]).sum() / 49.0
    return T


def ring_signal(T, qs=range(3, 9)):
    evals = np.linalg.eigvalsh((T + T.T) / 2)
    lam1 = evals[-1]
    r = abs(evals[-2] / lam1)
    # overflow-safe: F(q) = -q log lam1 - log sum_i (lam_i/lam1)^q
    ratios = evals / lam1
    Fs = np.array([-q * np.log(lam1) - np.log(np.sum(ratios ** q))
                   for q in qs])
    c1, c0 = np.polyfit(list(qs), Fs, 1)
    nonlin = float(np.abs(Fs - (c0 + c1 * np.array(list(qs)))).max())
    return r, c0, nonlin


def part_a():
    print("=" * 72)
    print("PART A -- required link correlation for O(1) curvature pricing")
    print("=" * 72)
    print(f"{'beta_eff':>9} {'r=|l2/l1|':>10} {'xi (links)':>11} "
          f"{'c0/edge':>9} {'nonlin':>9}  note")
    hits = {}
    for be in (1, 2, 4, 8, 16, 32, 64, 128):
        T = link_matrix(be)
        r, c0, nl = ring_signal(T)
        xi = -1.0 / np.log(r) if 0 < r < 1 else float("inf")
        sig = max(abs(c0), nl)
        note = ""
        for thr in (0.5, 2.2):
            if thr not in hits and sig >= thr:
                hits[thr] = (be, r, xi)
                note += f"  <- signal >= {thr}"
        print(f"{be:>9} {r:>10.4f} {xi:>11.2f} {c0:>+9.4f} {nl:>9.4f}{note}")
    print()
    for thr, (be, r, xi) in sorted(hits.items()):
        print(f"  REQUIREMENT (signal >= {thr}): beta_eff ~ {be}  "
              f"(x{be} amplification over the theory point), "
              f"r ~ {r:.3f}, xi ~ {xi:.1f} links")
    if not hits:
        print("  thresholds not reached in scanned range")
    return hits


# ---------------------------------------------------------------------------
# Part B: exact tilted heat-bath on a ring of q 7-state faces
# ---------------------------------------------------------------------------

def part_b(q=6, s_grid=(0.0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2), iters=800):
    print()
    print("=" * 72)
    print(f"PART B -- closure-maintaining history tilt, exact on a ring q={q}")
    print("=" * 72)
    T = link_matrix(1.0)
    v2 = -np.log(T)                      # pair potential defining V(x)
    shape = (7,) * q

    # V as a tensor: sum of pair potentials on the ring
    V = np.zeros(shape)
    for i in range(q):
        sl = [np.newaxis] * q
        sl[i] = slice(None)
        sl[(i + 1) % q] = slice(None)
        # build broadcastable pair term
        pair = v2
        expand = [np.newaxis] * q
        expand[i] = slice(None)
        expand[(i + 1) % q] = slice(None)
        # careful: axis order (i before or after (i+1)%q)
        if i < (i + 1) % q:
            V = V + pair[tuple(expand)]
        else:  # wrap: (q-1, 0) -> transpose
            V = V + pair.T[tuple(expand)]

    pi = np.exp(-V)
    pi /= pi.sum()
    sq = np.sqrt(pi)
    sq_safe = np.maximum(sq, 1e-300)

    def apply_K_sym(vec):
        """Symmetrized single-site heat-bath kernel, averaged over sites:
        H = D^{1/2} K D^{-1/2}, K_i f = E_pi[f | x_(not i)] (constant along
        axis i). Stationary eigenvector sqrt(pi) with eigenvalue 1."""
        f = vec / sq_safe
        out = np.zeros_like(vec)
        for i in range(q):
            cond = pi / pi.sum(axis=i, keepdims=True)
            avg = (cond * f).sum(axis=i, keepdims=True)
            out += sq * avg
        return out / q

    Vc = V - float((pi * V).sum())        # center the tilt observable
    results = []
    for s in s_grid:
        # tilted symmetric operator A = K_sym - s*Vc (+ shift for positivity)
        shift = s * float(Vc.max()) + 0.5
        vec = sq.copy()
        vec /= np.linalg.norm(vec)
        lam_prev = None
        for it in range(iters):
            nv = apply_K_sym(vec) - s * Vc * vec + shift * vec
            lam = float(np.linalg.norm(nv))
            vec = nv / lam
            if lam_prev is not None and abs(lam - lam_prev) < 1e-12 * lam:
                break
            lam_prev = lam
        mu = (vec * vec)
        mu /= mu.sum()
        # equal-time correlations of the label value m along the ring
        m_op = np.arange(-3, 4, dtype=float)
        grids = np.meshgrid(*([m_op] * q), indexing="ij")
        mean0 = (mu * grids[0]).sum()
        var0 = (mu * grids[0] ** 2).sum() - mean0 ** 2
        cs = []
        for d in range(1, q // 2 + 1):
            md = (mu * grids[0] * grids[d]).sum() \
                - mean0 * (mu * grids[d]).sum()
            cs.append(md / max(var0, 1e-12))
        # closure improvement the tilt buys, per face (its magnitude times s
        # lower-bounds the rarity rate of the selected histories)
        dV = float(((mu * V).sum() - (pi * V).sum())) / q
        results.append((s, cs, dV))
        print(f"  s={s:<5} corr(d)={['%+.3f' % c for c in cs]}  "
              f"dV/face = {dV:+.4f}  (rarity rate >= {abs(s * dV):.4f})")
    print("\n  corr(d): equal-time label correlation at ring distance d "
          "(untilted d=1 magnitude reflects r0~0.12);")
    print("  dV/face: closure improvement of the selected histories; "
          "s*|dV| lower-bounds their atypicality rate per face per unit "
          "refresh time.")
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=int, default=6)
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    part_a()
    part_b(q=args.q, iters=200 if args.fast else 800)
