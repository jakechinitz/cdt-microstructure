#!/usr/bin/env python3
"""FROM THEORY TOWARD REGGE, step 1: what effective geometric action does the
admissibility-closure weighting induce? (Analytic, transfer-matrix; no MC.)

THE QUESTION. 'Deriving Regge from the theory' means: expand the closure
term's geometry marginal -log E_labels[exp(-beta sum_tets E)] in geometric
quantities. Volume terms are known (the mu calibration). The first term that
can price CURVATURE comes from rings: the q tetrahedra around a slice edge,
where q IS the discrete curvature (deficit angle = 2*pi - q*theta). So the
target is the ring free energy F(q): its q-dependence is the induced
curvature-sector action.

THE CALCULATION. Consecutive tets in the ring share one face; out-of-ring
faces are integrated out uniformly (leading order of the cluster expansion;
corrections are O((lam2/lam1)^2) ~ 1%). The shared-face label is then a
7-state transfer chain:  T[m,m'] = (1/49) sum_{a,b} exp(-beta E(m,m',a,b)),
and F(q) = -log tr(T^q).  Decompose F(q) = c0 + c1*q + nonlinear:
  * c1 (per edge-tet incidence): sums to 6*c1*N_tet -- a VOLUME term,
    absorbed by the mu centering / k4. Not curvature.
  * c0 (per slice edge): prices N1 -- the curvature-LINEAR sector, i.e. the
    k0/Delta-type coupling the closure term induces on its own.
  * the nonlinear remainder: beyond-Regge (R^2-like); decays as
    (lam2/lam1)^q.

RESULT (beta=1, eta=eta*, lambda_inj=3 -- the theory point):

    lam2/lam1 = -0.120      (correlation range ~ half a lattice step)
    c0        = +0.019      per slice edge
    c1        = +0.528      per incidence (volume sector)
    nonlinear = ~1e-2 at q=3, ~5e-4 at q=5, negligible at physical q~5

  Against the bare couplings that sustain phase C (k0 ~ 2.2, Delta ~ 0.6):
  the closure weighting's induced curvature-sector coupling is ~1% of what
  is needed. With lambda_inj=0 (pure K^2, no injectivity) it drops another
  factor ~20 -- the injectivity constraint carries almost all of the
  geometric grip that exists.

WHAT THIS SETTLES.
  1. Closure-only ('moonshot', Regge removed) CANNOT produce 4D geometry
     from this weighting alone: the induced curvature pricing is two orders
     of magnitude too weak and the label correlations are too short-ranged
     to generate collective geometric order. This is now a computed number,
     not a judgment call.
  2. It PREDICTS the coupled sweeps should show only small geometric steers
     (percent-scale coupling shifts) -- a falsifiable cross-check against
     the stage-2 pilot: a large real-arm steer would CONTRADICT this
     analysis and demand explanation (artifact or beyond-cluster physics).
  3. For the emergence program, the load must be carried by the theory's
     OTHER dynamical ingredient -- the non-local memoryless refresh kernel
     (Many-Pasts), which the paper itself argues is required because local
     dynamics freeze (App. D.4). Coarse-graining the refresh dynamics, not
     the static weight, is the next math target.

Usage:  python induced_couplings.py [--beta 1.0 --eta ETA_STAR --lam 3.0]
"""
from __future__ import annotations
import argparse
import numpy as np

from v6_closure_run import build_energy_table, ETA_STAR


def ring_analysis(beta, eta, lam, qmax=8, verbose=True):
    E = build_energy_table(eta, lam)
    T = np.zeros((7, 7))
    for m in range(7):
        for mp in range(7):
            s = 0.0
            for a in range(7):
                for b in range(7):
                    s += np.exp(-beta * E[m, mp, a, b])
            T[m, mp] = s / 49.0
    evals = np.linalg.eigvalsh((T + T.T) / 2)
    lam1, lam2 = evals[-1], evals[-2]
    qs = np.arange(3, qmax + 1)
    Fs = np.array([-np.log(np.trace(np.linalg.matrix_power(T, q)))
                   for q in qs])
    c1, c0 = np.polyfit(qs, Fs, 1)
    nonlin = np.abs(Fs - (c0 + c1 * qs))
    if verbose:
        print(f"[beta={beta} eta={eta:.6g} lambda={lam}]")
        print(f"  spectrum: lam1={lam1:.6f} lam2={lam2:.6f} "
              f"ratio={lam2 / lam1:+.4f}")
        for q, F, r in zip(qs, Fs, nonlin):
            print(f"  q={q}: F={F:.4f}  nonlinear={r:.2e}")
        print(f"  induced per-edge (curvature-sector) coupling c0 = {c0:+.4f}")
        print(f"  induced per-incidence (volume-sector)   c1 = {c1:+.4f}")
        print(f"  bare couplings sustaining phase C: k0 ~ 2.2, Delta ~ 0.6"
              f"  ->  c0 is ~{abs(c0) / 2.2 * 100:.1f}% of k0")
    return {"lam_ratio": lam2 / lam1, "c0": c0, "c1": c1,
            "nonlin_max": float(nonlin.max())}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=ETA_STAR)
    ap.add_argument("--lam", type=float, default=3.0)
    args = ap.parse_args()
    ring_analysis(args.beta, args.eta, args.lam)
    print()
    ring_analysis(args.beta, args.eta, 0.0)     # no-injectivity comparison
