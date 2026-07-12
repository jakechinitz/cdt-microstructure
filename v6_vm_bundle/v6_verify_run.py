#!/usr/bin/env python3
"""VERIFICATION / known-answer run: bare-Regge 4D CDT on the v6 engine.

PURPOSE: prove the v6 engine reproduces KNOWN 4D CDT physics before any theory
is added. We run pure Einstein-Regge dynamics at the canonical de Sitter point
and check for the literature signatures of phase C:

  * Hausdorff dimension d_H climbing toward ~4 (vs torus rails: a true 4-torus
    of the same size reads d_H ~ 3.1 with this estimator; 2-torus ~1.9).
  * Spatial-volume profile V3(tau) forming a localized "blob" (de Sitter
    cos^3 universe) instead of a flat minimal "stalk".

If this passes, the engine is validated AND it is itself a demonstration that
"the graph generates 4D spacetime". Only then is the with-theory run
(v6_theory_run.py) meaningful.

CANONICAL SETTINGS (Ambjorn-Gorlich-Jurkiewicz-Loll; Gorlich PhD thesis):
  (kappa_0, Delta) = (2.2, 0.6), kappa_4 ~ 0.9 (pseudo-critical), fix N_41,
  S^1 x S^3, T ~ 80 time slices. De Sitter is a LARGE-volume phenomenon: the
  literature needs N_41 >~ 20k-160k (N_4 ~ 45k-360k); below ~ a few x10^4 the
  universe collapses to the stalk. So target N_41 as high as the machine allows.

USAGE (drop on any VM with numpy/scipy):
  python v6_verify_run.py --target-n41 40000 --K 64 --checkpoint ckpt.json
  python v6_verify_run.py --resume ckpt.json   # continue a long run
"""
from __future__ import annotations
import argparse
from v6_run_lib import (run, dual_adjacency, hausdorff_dim, volume_profile,
                        profile_metrics, torus_rails)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target-n41", type=int, default=40000,
                   help="target N_41 (de Sitter needs >~ 20000-160000)")
    p.add_argument("--K", type=int, default=64, help="number of time slices (T~80 canonical)")
    p.add_argument("--k0", type=float, default=2.2)
    p.add_argument("--Delta", type=float, default=0.6)
    p.add_argument("--k4", type=float, default=0.9)
    p.add_argument("--eps", type=float, default=1e-4, help="N_41 volume-fix strength (small!)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-sweeps", type=int, default=200000)
    p.add_argument("--measure-every", type=int, default=50)
    p.add_argument("--checkpoint", type=str, default="v6_verify_ckpt.json")
    p.add_argument("--resume", type=str, default=None)
    p.add_argument("--wall-hours", type=float, default=None, help="stop after this many hours")
    p.add_argument("--causal-slices", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="enforce the CDT foliation (slices stay closed 3-manifolds; "
                        "standard AJL ensemble). --no-causal-slices reproduces the "
                        "pre-fix generalized ensemble.")
    p.add_argument("--tune-k4", type=int, default=0, metavar="ROUNDS",
                   help="auto-tune k4 to pseudo-criticality over N short "
                        "bursts before the main run (phase_grid.py protocol: "
                        "the eps-pin residual is the feedback error, "
                        "k4 += 2*eps*(N41-target)). Fixed k4=0.9 was found "
                        "to fight the pin across the grid; production runs "
                        "at the tuned point should use --tune-k4 3.")
    p.add_argument("--tune-burst", type=int, default=60,
                   help="sweeps per tuning burst (after the first burst, "
                        "which also grows to the target volume)")
    args = p.parse_args()

    # --- optional k4 auto-tuning (pseudo-criticality via the pin residual) --
    k4 = args.k4
    resume = args.resume
    if args.tune_k4 > 0:
        for rd in range(1, args.tune_k4 + 1):
            T = run(f"k4-tune round {rd}/{args.tune_k4}", k0=args.k0,
                    Delta=args.Delta, k4=k4, target_N41=args.target_n41,
                    K=args.K, eps=args.eps, seed=args.seed,
                    max_sweeps=args.tune_burst, measure_every=args.tune_burst,
                    checkpoint=args.checkpoint, resume=resume,
                    causal=args.causal_slices, verbose=False)
            resume = args.checkpoint
            n41 = T.type_counts()[0]
            dk = 2 * args.eps * (n41 - args.target_n41)
            k4 += dk
            print(f"# [k4-tune] round {rd}: N41={n41} (target "
                  f"{args.target_n41})  pin residual => dk4={dk:+.4f}  "
                  f"k4 -> {k4:.4f}", flush=True)
        print(f"# [k4-tune] final k4 = {k4:.4f}", flush=True)

    T = run("bare-Regge (verification)", k0=args.k0, Delta=args.Delta, k4=k4,
            target_N41=args.target_n41, K=args.K, eps=args.eps, seed=args.seed,
            max_sweeps=args.max_sweeps, measure_every=args.measure_every,
            checkpoint=args.checkpoint, resume=resume,
            causal=args.causal_slices,
            wall_budget_s=(args.wall_hours * 3600 if args.wall_hours else None))

    # final verdict
    ids, adj = dual_adjacency(T)
    dH = hausdorff_dim(adj)
    prof = volume_profile(T)
    pm = profile_metrics(prof)
    rails = torus_rails(T.n_pent())
    okf, repf = getattr(T, "_final_verify", (None, {}))
    ga = repf.get("gluing", {})
    print("\n" + "=" * 64)
    print("  VERIFICATION VERDICT")
    print("=" * 64)
    print(f"  final N4={T.n_pent()}  N41={T.type_counts()[0]}")
    print(f"  manifold check (gluing-based + S^3 links): "
          f"{'PASS' if okf else 'FAIL'}  "
          f"[gluing_ok={repf.get('ok_gluing_only')}, "
          f"links={repf.get('link_failures')}, "
          f"simplicial={ga.get('is_simplicial')}, chi={repf.get('euler_char')}]")
    print(f"  d_H = {dH:.2f}")
    print(f"  rails (matched size): 2-torus={rails[2][1]:.2f}  "
          f"3-torus={rails[3][1]:.2f}  4-torus={rails[4][1]:.2f}")
    print("  --- de Sitter volume-profile diagnostics ---")
    print(f"  blob score    = {pm['blob_score']:.2f}   (1=stalk/flat, >>1=blob)")
    print(f"  active slices = {pm['active_slices']} / {T.K}   "
          f"(blob extent; full stalk-collapse ~ 1-2)")
    print(f"  max slice     = {pm['max_slice']:.0f}   mean stalk = {pm['mean_stalk']:.1f}")
    if pm['cos3_relerr'] is not None:
        print(f"  cos^3 fit     : width={pm['cos3_width']:.1f} slices  "
              f"rel.RMS err={pm['cos3_relerr']:.3f}  (lower => more de-Sitter-like)")
    print(f"  centered V3(tau) = {pm['centered']}")
    near4 = abs(dH - rails[4][1]) < abs(dH - rails[2][1])
    blob = pm['blob_score'] > 1.5 and pm['active_slices'] >= 3
    print(f"\n  VERDICT: d_H is closer to the "
          f"{'4-TORUS (4D emergent!)' if near4 else '2-torus (still collapsed -- grow larger)'}")
    print(f"           volume profile {'shows a localized blob (de Sitter-like)' if blob else 'is flat/stalk (collapsed -- grow larger)'}")


if __name__ == "__main__":
    main()
