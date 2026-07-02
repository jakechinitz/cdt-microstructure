#!/usr/bin/env python3
"""EPRL term diagnostic: is the (centered) vertex-amplitude term INERT or does
it have TEETH?

The matched-volume beta sweep can only be interpreted if the centered EPRL term
actually carries variance across configurations. If the centered
per-pentachoron log-amplitude is a near-delta spike, the term is inert: every
"flat d_H vs beta" result afterward is an artifact of a dead term, not evidence
the theory tolerates 4D. If it has real spread, a flat result is a genuine
finding ("the term fluctuates but does not steer dimension").

This is a ZERO-SWEEP test: it builds one thermalized-label config and histograms
  centered_p = -log|A_v(faces of p)| - mu ,   mu = mean over pentachora.

Reported:
  * sigma = std of centered contributions (inert if <0.05, has-teeth if >0.5)
  * the effective per-move action scale beta*sigma at the swept betas
  * whether the variance aligns with the (4,1)/(3,2) simplex-type split -- the
    geometric DOF that actually controls the volume profile / d_H. Variance that
    is large but type-blind means the term has teeth largely ORTHOGONAL to the
    geometry, which is itself the correct reading of a flat sweep.

Usage:
  python eprl_term_diagnostic.py [--K 16 --grow 8000 --hb 15 --seed 0]
"""
from __future__ import annotations
import argparse
import numpy as np

from v6_cdt import build_s1xs3
from v6_cdt_run import propose_and_apply
from v6_theory_run import Intertwiners, make_heatbath
from vertex_tensor import FaithfulVertex


def per_pent_neglog(T, intw, Ttensor):
    out, typ = [], []
    for p in T.pent:
        faces = tuple(intw.lab[frozenset((p, T.nbr[p][i]))] for i in range(5))
        out.append(-np.log(max(abs(Ttensor[faces]), 1e-300)))
        a, b = T.ptype(p)
        typ.append("41" if {a, b} == {4, 1} else "32" if {a, b} == {3, 2} else "x")
    return np.array(out), np.array(typ)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vertex", default="vertex_j3.npz")
    ap.add_argument("--K", type=int, default=16)
    ap.add_argument("--grow", type=int, default=8000, help="proposal attempts to grow/thermalize geometry")
    ap.add_argument("--hb", type=int, default=15, help="heat-bath passes to equilibrate labels")
    ap.add_argument("--hb-beta", type=float, default=1.0,
                    help="label exponent for the heat-bath. Use the run's "
                         "beta_eprl to see the term at that coupling (the "
                         "production heat-bath is now beta-consistent).")
    ap.add_argument("--symmetrize", action="store_true",
                    help="slot-symmetrize the tensor first (matches the "
                         "production default in v6_theory_run.py)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--betas", default="0.05 0.1 0.3")
    args = ap.parse_args()

    Tt = FaithfulVertex.load(args.vertex).dense_tensor().astype(np.float64)
    if args.symmetrize:
        from v6_theory_run import slot_symmetrize
        Tt = slot_symmetrize(Tt)
        print("[tensor] slot-symmetrized (production default)")
    D = Tt.shape[0]
    rng = np.random.default_rng(args.seed)

    allraw = -np.log(np.abs(Tt).clip(1e-300))
    print(f"[tensor] D={D}  -log|amp| over all {D**5} entries: "
          f"min={allraw.min():.2f} max={allraw.max():.2f} "
          f"mean={allraw.mean():.2f} std={allraw.std():.2f}")

    T = build_s1xs3(K=args.K)
    for _ in range(args.grow):
        propose_and_apply(T, rng, 2)
    intw = Intertwiners(T, D, rng)
    hb = make_heatbath(intw, Tt, beta=args.hb_beta)
    for _ in range(args.hb):
        hb(T, rng)

    raw, typ = per_pent_neglog(T, intw, Tt)
    mu = raw.mean()
    cen = raw - mu
    sigma = cen.std()
    print(f"\n[heat-bath equilibrated at exponent {args.hb_beta}, N4={len(raw)}]"
          f"  mean cost = {mu:.4f}  (NB: the run's centering mu is the free "
          f"energy from calibrate_mu_ti, not this mean)")
    print(f"  centered -log|amp|: sigma={sigma:.4f}  min={cen.min():.3f} max={cen.max():.3f}")
    qs = np.percentile(cen, [1, 10, 25, 50, 75, 90, 99])
    print(f"  percentiles [1,10,25,50,75,90,99] = {np.round(qs, 3)}")

    intw_u = Intertwiners(T, D, rng)
    raw_u, _ = per_pent_neglog(T, intw_u, Tt)
    print(f"  (uniform-random labels sigma = {(raw_u-raw_u.mean()).std():.4f}  "
          f"-> heat-bath narrows the distribution)")

    m41, m32 = cen[typ == "41"], cen[typ == "32"]
    if len(m41) and len(m32):
        split = abs(m41.mean() - m32.mean())
        print(f"\n[geometry coupling] mean centered by simplex type: "
              f"(4,1)={m41.mean():+.4f}  (3,2)={m32.mean():+.4f}  |split|={split:.4f}")
        print(f"  -> {'type-BLIND: teeth largely orthogonal to geometry (flat sweep = genuine)' if split < 0.2*sigma else 'type-COUPLED: term distinguishes the geometric DOF'}")

    print(f"\n[verdict] sigma={sigma:.4f}  "
          f"{'INERT (dead term)' if sigma < 0.05 else 'HAS TEETH' if sigma > 0.5 else 'WEAK'}")
    for b in (float(x) for x in args.betas.split()):
        print(f"  beta={b}: per-pent action scale beta*sigma = {b*sigma:.4f}")


if __name__ == "__main__":
    main()
