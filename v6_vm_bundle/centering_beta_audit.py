#!/usr/bin/env python3
"""Centering audit: is the 'volume-neutral' EPRL term actually volume-neutral?

The v6 theory run subtracts a fixed per-pentachoron constant mu from the raw
EPRL cost -log|A_v| so the coupled term measures fluctuations only. mu is
auto-calibrated ONCE, on the FIRST action evaluation -- which, on a fresh run
(or a resume from a bare checkpoint), happens with UNIFORM-RANDOM intertwiner
labels, i.e. mu ~= mean(-log|A|) over all tensor entries (~30.53 for
vertex_j3.npz). The per-sweep heat-bath then concentrates labels on
high-amplitude entries, so the equilibrated mean cost is LOWER than mu, and the
"centered" term becomes an EXTENSIVE offset

    S_EPRL / N4  ~=  (equilibrated mean cost) - mu  <  0 ,

i.e. an uncontrolled effective cosmological-constant shift

    k4_eff = k4 + beta * offset(beta)

that grows monotonically with beta. Since the volume penalty pins N41 only,
the surplus drives the UNPINNED N32 sector, so a beta sweep at 'matched volume'
doubles as an uncontrolled k4/N32 sweep -- a size->shape confound in exactly
the direction the sweep is trying to measure.

A second, compounding defect: the shipped heat-bath (make_heatbath /
v6_theory_run.py) resamples labels with weight prod|A|^1 REGARDLESS of beta,
while the geometry moves feel beta * S_EPRL. The joint chain therefore has no
consistent stationary measure, and the labels always sit at the beta=1
(maximally concentrated) conditional -- which makes the offset above the FULL
-2.5 per pentachoron at every beta>0 instead of the much smaller
beta-consistent value.

This script measures offset(beta) with a beta-CONSISTENT heat-bath on a fixed
small geometry, and prints the resulting effective-k4 shifts for both the
correct (beta-consistent) and the shipped (beta=1) label distributions.

Reference output (vertex_j3.npz, K=16, grow=8000, 12 passes, seed 0/1):

     beta  mean cost    std   offset(beta)   k4 shift (correct)  k4 shift (shipped, beta=1 labels)
      0.0     30.68   1.87      +0.14            +0.000              -0.000
     0.05     30.43   1.79      -0.10            -0.005              -0.128
      0.1     30.22   1.62      -0.32            -0.032              -0.256
      0.3     29.71   1.38      -0.82            -0.247              -0.769
      1.0     27.97   0.89      -2.56            -2.563              -2.563

Usage:  python centering_beta_audit.py [--K 16 --grow 8000 --hb 12 --seed 0]
"""
from __future__ import annotations
import argparse
import numpy as np

from v6_cdt import build_s1xs3
from v6_cdt_run import propose_and_apply
from v6_theory_run import Intertwiners


def heatbath_beta(T, intw, Tt, beta, passes, rng, D):
    """Same conditional as make_heatbath() but at the CORRECT exponent beta."""
    for _ in range(passes):
        for key in list(intw.lab.keys()):
            ps = [p for p in key if p in T.pent]
            if len(ps) != 2:
                continue
            logw = np.zeros(D)
            for c in range(D):
                intw.lab[key] = c
                lw = 0.0
                for p in ps:
                    faces = tuple(intw.lab[frozenset((p, T.nbr[p][i]))]
                                  for i in range(5))
                    lw += np.log(Tt[faces])
                logw[c] = beta * lw
            logw -= logw.max()
            w = np.exp(logw)
            w /= w.sum()
            intw.lab[key] = int(rng.choice(D, p=w))


def mean_cost(T, intw, Tt):
    cs = []
    for p in T.pent:
        faces = tuple(intw.lab[frozenset((p, T.nbr[p][i]))] for i in range(5))
        cs.append(-np.log(Tt[faces]))
    return float(np.mean(cs)), float(np.std(cs))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vertex", default="vertex_j3.npz")
    ap.add_argument("--K", type=int, default=16)
    ap.add_argument("--grow", type=int, default=8000)
    ap.add_argument("--hb", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--betas", default="0.0 0.05 0.1 0.3 1.0")
    args = ap.parse_args()

    Tt = np.abs(np.asarray(np.load(args.vertex)["tensor"], float)).clip(1e-300)
    D = Tt.shape[0]
    mu_uniform = float(-np.log(Tt).mean())
    print(f"[tensor] uniform-label mean cost (== the auto-calibrated mu on a "
          f"fresh/bare-resumed run): {mu_uniform:.3f}")

    rng = np.random.default_rng(args.seed)
    T = build_s1xs3(K=args.K)
    for _ in range(args.grow):
        propose_and_apply(T, rng, 2)
    print(f"[geometry] N4={T.n_pent()}  (fixed for all betas)\n")

    betas = [float(x) for x in args.betas.split()]
    # shipped behaviour = labels always at the beta=1 conditional
    intw1 = Intertwiners(T, D, np.random.default_rng(args.seed + 1))
    heatbath_beta(T, intw1, Tt, 1.0, args.hb, rng, D)
    m1, _ = mean_cost(T, intw1, Tt)
    off1 = m1 - mu_uniform

    print(f"{'beta':>6} {'mean cost':>10} {'std':>6} {'offset':>8} "
          f"{'k4 shift (correct)':>19} {'k4 shift (shipped)':>19}")
    for beta in betas:
        intw = Intertwiners(T, D, np.random.default_rng(args.seed + 1))
        heatbath_beta(T, intw, Tt, beta, args.hb, rng, D)
        m, s = mean_cost(T, intw, Tt)
        off = m - mu_uniform
        print(f"{beta:>6} {m:>10.3f} {s:>6.3f} {off:>+8.3f} "
              f"{beta * off:>+19.4f} {beta * off1:>+19.4f}")
    print(f"\n(bare k4 in the sweep = 0.9; any |k4 shift| that is a sizable "
          f"fraction of that is an uncontrolled volume/N32 drive, not theory.)")


if __name__ == "__main__":
    main()
