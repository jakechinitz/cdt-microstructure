#!/usr/bin/env python3
"""Criterion 3a: the per-move EPRL action contribution at beta -- the diagnostic
that decides GOOD null vs EMPTY null (see NULL_INTERPRETATION_ARGUMENT.md).

A flat blob is only the GOOD (publishable) null if the EPRL term actually entered
accept/reject decisions -- i.e. its per-move action change beta*dS_EPRL was O(1)
(comparable to the Regge per-move change) for a non-negligible fraction of moves.
If beta*dS_EPRL was sub-threshold on essentially every move, the term washed out
and a flat result is the EMPTY null (a statement about our setup, not gravity).

This measures, on a thermalized config, the distribution of |beta * dS_EPRL| over
attempted Pachner moves, and compares it to |dS_Regge| on the same moves. Reports
the fraction of moves where the EPRL term was dynamically relevant (>~0.1 and
>~1.0 on the Metropolis scale, and relative to Regge).

Usage (from v6_vm_bundle/, after thermalization):
  .venv/bin/python3 per_move_eprl_scale.py --beta 0.3 --K 24 --grow 40000 --moves 4000
"""
from __future__ import annotations
import argparse
import numpy as np

from v6_cdt import build_s1xs3, regge_action
from v6_cdt_run import propose_and_apply
from v6_theory_run import Intertwiners, make_heatbath, IncrementalEPRL, Centering
from vertex_tensor import FaithfulVertex


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vertex", default="vertex_j3.npz")
    ap.add_argument("--beta", type=float, default=0.3)
    ap.add_argument("--K", type=int, default=24)
    ap.add_argument("--grow", type=int, default=40000)
    ap.add_argument("--hb", type=int, default=12)
    ap.add_argument("--moves", type=int, default=4000, help="attempted moves to sample dS over")
    ap.add_argument("--k0", type=float, default=2.2)
    ap.add_argument("--Delta", type=float, default=0.6)
    ap.add_argument("--k4", type=float, default=0.9)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    Tt = FaithfulVertex.load(args.vertex).dense_tensor().astype(np.float64)
    D = Tt.shape[0]
    rng = np.random.default_rng(args.seed)
    T = build_s1xs3(K=args.K)
    for _ in range(args.grow):
        propose_and_apply(T, rng, 2)
    intw = Intertwiners(T, D, rng)
    hb = make_heatbath(intw, Tt)
    for _ in range(args.hb):
        hb(T, rng)

    centering = Centering(enabled=True)
    inc = IncrementalEPRL(intw, Tt, "regge_plus_eprl", args.k0, args.Delta,
                          args.k4, args.beta, rng, centering=centering)
    inc.full(T)

    eprl_ds, regge_ds = [], []
    for _ in range(args.moves):
        regge_before = regge_action(T, args.k0, args.Delta, args.k4)
        T.begin_record()
        mt, ok, undo = propose_and_apply(T, rng, 2)
        if not ok:
            T.take_record()
            continue
        added, removed = T.take_record()
        dS_total = inc.delta(T, added, removed)          # beta-weighted full geom delta
        regge_after = regge_action(T, args.k0, args.Delta, args.k4)
        dS_regge = regge_after - regge_before
        dS_eprl = dS_total - dS_regge                    # the beta*S_EPRL part
        eprl_ds.append(abs(dS_eprl))
        regge_ds.append(abs(dS_regge))
        # reject-and-restore so we keep sampling near the thermalized state
        T.begin_record(); undo(); ua, ur = T.take_record(); inc.reject(T, ua, ur)

    e = np.array(eprl_ds); r = np.array(regge_ds)
    print(f"[per-move EPRL scale]  beta={args.beta}  N4={T.n_pent()}  "
          f"sampled {len(e)} accepted-then-undone moves")
    print(f"  |beta*dS_EPRL|: median={np.median(e):.3f}  mean={e.mean():.3f}  "
          f"90th={np.percentile(e,90):.3f}  max={e.max():.3f}")
    print(f"  |dS_Regge|    : median={np.median(r):.3f}  mean={r.mean():.3f}")
    frac_01 = np.mean(e > 0.1)
    frac_1 = np.mean(e > 1.0)
    frac_vs_regge = np.mean(e > 0.5 * np.maximum(r, 1e-9))
    print(f"  fraction of moves with |beta*dS_EPRL| > 0.1 (Metropolis-relevant) = {frac_01:.2f}")
    print(f"  fraction with > 1.0 (strongly relevant)                          = {frac_1:.2f}")
    print(f"  fraction where EPRL term >= 0.5*|dS_Regge| (competes w/ gravity)  = {frac_vs_regge:.2f}")
    print()
    # The verdict is NOT binary. Three regimes, and "loud" splits two ways:
    if frac_01 < 0.2:
        print("  REGIME: WASHED OUT -- EPRL term sub-threshold on ~all moves. A flat blob")
        print("    is the EMPTY null by weakness; fix is higher beta.")
    else:
        print("  REGIME: term is DYNAMICALLY PRESENT (not washed out) -- it enters")
        print("    accept/reject. This rules out empty-null-by-weakness. BUT 'present' is")
        print("    not yet 'GOOD null': two further checks REQUIRED, not handled here --")
        print(f"    (i) ORTHOGONALITY: with type-split/sigma~0.06, a loud term may still be")
        print(f"        geometry-blind (loud NOISE, not steering). A flat blob is the GOOD")
        print(f"        null only if you argue the r_nn=0.14 aligned channel had real")
        print(f"        leverage -- see NULL_INTERPRETATION_ARGUMENT.md.")
        # crude acceptance estimate from the sampled dS (Metropolis cap at 1)
        acc = float(np.mean(np.minimum(1.0, np.exp(-e))))
        print(f"    (ii) FREEZING: a loud term suppresses acceptance. Rough accept on these")
        print(f"        moves ~ <min(1,exp(-|beta*dS_EPRL|))> = {acc:.3f}. If this is tiny,")
        print(f"        beta=0.3 may be a FROZEN chain (prereg gate G5) -- itself not a")
        print(f"        result. Confirm the live beta=0.3 run is still moving (N4/blob")
        print(f"        evolving) before trusting any flat read.")


if __name__ == "__main__":
    main()
