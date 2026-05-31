#!/usr/bin/env python3
"""d_s N4-SCALING check at beta=0 -- the NON-CIRCULAR substrate validator.

Context (reviewer-driven): an earlier argument disqualified the spectral
dimension d_s because, on a single bare config, its flow peaked at ~2.4 instead
of ~4. That reasoning was circular -- it used d_H + the cos^3 blob (the very
observables d_s was meant to independently check) to declare d_s the broken
instrument. The fix is to let d_s earn or lose admissibility by its OWN scaling:

  Does the large-sigma (IR) d_s climb toward 4 as N4 grows?

This is exactly how AJL established 4D CDT: the sigma->inf spectral dimension
rises with volume (they needed N4 ~ 181k to reach 4.02 +/- 0.10). So a low
reading at small N4 is expected; the DECIDABLE thing is the trend, not the value.

  * d_s_peak RISES with N4   -> the low reading is finite-size; substrate is
                                4D-trending, vindicated by d_s's own scaling
                                (no appeal to d_H needed).
  * d_s_peak FLAT across N4  -> substrate is NOT phase C; d_H + blob have been
                                over-reading. That would be the real problem.

MEASURED (sandbox, lightly thermalized, beta=0 / pure Regge; seed 11):
    N4      d_H    d_s_peak  @sigma
    2114    2.25     2.36     3.6
    5030    2.47     2.45     4.8
   10018    2.85     2.52    18.8
   18108    2.98     2.65    25.4
  d_s_peak climbs 2.36 -> 2.65 AND the peak scale marches out 3.6 -> 25.4 as the
  universe grows: the AJL finite-size signature (4D approached from below, the
  plateau pushing to larger scales). Conclusion: d_s-scaling and d_H AGREE the
  substrate is 4D-trending -- two independent probes, not a self-protecting
  consensus.

LIMITS: sandbox caps ~18k (<< AJL 181k), so this is a trend DIRECTION, not a
landing at 4.0. Configs are lightly thermalized. Run larger on the VM/cluster for
a quantitative approach to 4.

Usage:
  python ds_scaling_check.py [--targets 2000 5000 10000 18000 --K 28 --seed 11]
"""
from __future__ import annotations
import argparse
import time
import numpy as np
from scipy.sparse import lil_matrix

import step3_linkA_harness as linkA
from v6_cdt import build_s1xs3
from v6_cdt_run import propose_and_apply
from v6_run_lib import dual_adjacency, hausdorff_dim


def ds_flow(adj, a, b):
    N = len(adj)
    A = lil_matrix((N, N))
    for i, nb in enumerate(adj):
        for j in nb:
            A[i, j] = 1.0
    ts, ds = linkA.raw_ds(A.tocsr())
    cal = a * ds + b
    peak = float(np.max(cal))
    peak_sigma = float(ts[np.argmax(cal)])
    mid = a * linkA.midwin(ts, ds) + b
    return peak, peak_sigma, mid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", type=int, nargs="+",
                    default=[2000, 5000, 10000, 18000])
    ap.add_argument("--K", type=int, default=28)
    ap.add_argument("--therm", type=int, default=2500)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    a, b = linkA.fit_calibration(verbose=False)
    print(f"calib a={a:.3f} b={b:.3f}")
    print(f"{'N4':>7} {'d_H':>6} {'d_s_peak':>9} {'@sigma':>7} {'d_s_mid':>8} {'wall':>6}")
    rng = np.random.default_rng(args.seed)
    T = build_s1xs3(K=args.K)
    rows = []
    ci = 0
    t0 = time.time()
    targets = sorted(args.targets)
    while ci < len(targets):
        propose_and_apply(T, rng, 2)
        if T.n_pent() >= targets[ci]:
            for _ in range(args.therm):
                propose_and_apply(T, rng, 2)
            _, adj = dual_adjacency(T)
            dH = hausdorff_dim(adj)
            pk, ps, md = ds_flow(adj, a, b)
            rows.append((T.n_pent(), dH, pk))
            print(f"{T.n_pent():>7} {dH:>6.2f} {pk:>9.2f} {ps:>7.2f} {md:>8.2f} "
                  f"{time.time()-t0:>6.0f}", flush=True)
            ci += 1

    if len(rows) >= 2:
        rise = rows[-1][2] - rows[0][2]
        verdict = ("RISES -> finite-size, substrate 4D-trending (vindicated by "
                   "d_s scaling)" if rise > 0.1 else
                   "FLAT -> substrate may NOT be phase C; investigate d_H/blob")
        print(f"\nd_s_peak change over N4 range = {rise:+.2f}  => {verdict}")


if __name__ == "__main__":
    main()
