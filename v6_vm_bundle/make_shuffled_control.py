#!/usr/bin/env python3
"""Make an entry-SHUFFLED copy of the vertex tensor: the placebo control.

The shuffled tensor has exactly the same 16807-entry value distribution as
vertex_j3.npz (same mean/std of -log|A|, same max|A|, same 'teeth') but ALL
EPRL structure -- every correlation between intertwiner labels -- is destroyed.

Why this control matters: the theory sweep asks whether the EPRL amplitude
steers the CDT geometry. Every artifact identified in the coupling machinery
(centering offset, beta-less heat-bath, uncorrected label births, arbitrary
slot assignment) depends only on the ENTRY STATISTICS of the tensor, not on
its EPRL structure. So:

  * If a beta>0 run with the SHUFFLED tensor reproduces the same d_H / blob /
    N32 movement as the real tensor at the same beta, the observed steer is a
    machinery artifact and says nothing about EPRL.
  * Only a DIFFERENCE between real-tensor and shuffled-tensor runs at the same
    beta can be attributed to the EPRL amplitude itself.

First small-volume result (K=16, N41=2000, eps=0.01, k4=0.9, 400 sweeps,
matched N41, this repo, 2026-07): the shuffled run REPRODUCED the real-tensor
beta=0.3 signature almost exactly (final N4 2660 vs real-tensor 2714 vs
control 2962; final blob 3.20 vs 3.21 vs 2.70) -- i.e. at this size the
entire visible 'steer' is machinery, not EPRL. Rerun at production volume
before quoting. Details: SIM_AUDIT_coupling_misspecifications.md.

Usage:
  python make_shuffled_control.py [--in vertex_j3.npz] [--out vertex_j3_shuffled.npz] [--seed 12345]
  python v6_theory_run.py --vertex vertex_j3_shuffled.npz --beta-eprl 0.3 ...
"""
from __future__ import annotations
import argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", default="vertex_j3.npz")
    ap.add_argument("--out", dest="dst", default="vertex_j3_shuffled.npz")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    d = dict(np.load(args.src))
    T = np.asarray(d["tensor"], float)
    flat = T.ravel().copy()
    np.random.default_rng(args.seed).shuffle(flat)
    d["tensor"] = flat.reshape(T.shape)
    np.savez(args.dst, **d)
    print(f"[shuffled control] {args.src} -> {args.dst}  (seed={args.seed})")
    print(f"  max|A| preserved: {np.abs(d['tensor']).max():.6e}")
    print(f"  NOTE: metadata (j, validated, ...) copied verbatim so the loader "
          f"accepts it; this file is a PLACEBO, never a physics tensor.")


if __name__ == "__main__":
    main()
