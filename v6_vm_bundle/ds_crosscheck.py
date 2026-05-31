#!/usr/bin/env python3
"""Spectral-dimension (d_s) cross-check for v6 CDT dual graphs -- EXPLICITLY
NOT VALIDATED as a standalone verdict instrument. Read this header before use.

The v6 engine deliberately uses d_H (Hausdorff, shell-growth) + the volume
profile, NOT d_s, because the diffusion d_s estimator is known to misread CDT
dual graphs (an extended 3D CDT reads d_s ~ 2). This module wires the validated-
on-tori linkA heat-kernel estimator (step3_linkA_harness) onto a v6 dual graph
so d_s can be reported ALONGSIDE d_H as a qualitative cross-check.

WHY IT IS NOT A VERDICT INSTRUMENT (measured evidence, sandbox volume scan):
    N4      d_H (trusted)   d_s_cal (linkA, tori-calibrated)
    1434       1.95              2.22
    2202       2.31              2.26
    8038       2.77              2.48   <-- d_s now LAGS d_H
  d_s and d_H agree at small N but DIVERGE as volume grows: d_s under-reads,
  exactly as the v6 warning predicts. So d_s here tracks QUALITATIVELY (both rise
  with volume) but cannot be trusted QUANTITATIVELY. Do not gate a THEORY_LIVES /
  THEORY_DIES verdict on d_s. The trustworthy second observable is the volume
  profile / cos^3 blob (already in v6_run_lib.profile_metrics).

Legitimate uses:
  * report d_s next to d_H as a consistency cross-check ("do both move the same
    way across beta?"), with the divergence caveat stated;
  * run validate_against_dH() to re-measure the d_s-vs-d_H relationship at the
    volume you actually care about before quoting anything.

Usage:
  from ds_crosscheck import ds_of_dual, validate_against_dH
  d_s = ds_of_dual(T)                      # calibrated d_s of a triangulation
  validate_against_dH(targets=[...])       # re-run the divergence check
"""
from __future__ import annotations
import numpy as np
from scipy.sparse import lil_matrix

import step3_linkA_harness as linkA
from v6_run_lib import dual_adjacency, hausdorff_dim


def _sparse_adj(adj):
    N = len(adj)
    A = lil_matrix((N, N))
    for i, nb in enumerate(adj):
        for j in nb:
            A[i, j] = 1.0
    return A.tocsr()


# Cache the tori calibration (a,b); it is geometry-independent.
_CAL = None


def calibration():
    global _CAL
    if _CAL is None:
        _CAL = linkA.fit_calibration(verbose=False)
    return _CAL


def ds_of_dual(T, a=None, b=None):
    """Calibrated linkA d_s of triangulation T's dual graph.
    NOT a verdict number -- see module header (under-reads vs d_H at large N)."""
    if a is None or b is None:
        a, b = calibration()
    _, adj = dual_adjacency(T)
    ts, ds = linkA.raw_ds(_sparse_adj(adj))
    return a * linkA.midwin(ts, ds) + b


def validate_against_dH(targets=(2000, 8000, 20000, 45000), K=24, seed=3,
                        therm=3000, verbose=True):
    """Re-measure d_s vs d_H across a volume scan on freshly grown CDT dual
    graphs. Prints the pair at each target so you can judge, at YOUR volume,
    whether d_s tracks d_H closely enough to quote. Returns list of dicts."""
    from v6_cdt import build_s1xs3
    from v6_cdt_run import propose_and_apply
    a, b = calibration()
    rng = np.random.default_rng(seed)
    T = build_s1xs3(K=K)
    out = []
    ci = 0
    targets = list(targets)
    if verbose:
        print(f"calibration a={a:.3f} b={b:.3f}")
        print(f"{'N4':>7} {'d_H':>6} {'d_s_cal':>8} {'gap':>6}")
    for _ in range(10_000_000):
        propose_and_apply(T, rng, 2)
        if ci < len(targets) and T.n_pent() >= targets[ci]:
            for _ in range(therm):
                propose_and_apply(T, rng, 2)
            _, adj = dual_adjacency(T)
            dH = hausdorff_dim(adj)
            ds = ds_of_dual(T, a, b)
            rec = {"N4": T.n_pent(), "d_H": dH, "d_s": ds, "gap": ds - dH}
            out.append(rec)
            if verbose:
                print(f"{rec['N4']:>7} {dH:>6.2f} {ds:>8.2f} {ds-dH:>+6.2f}", flush=True)
            ci += 1
            if ci >= len(targets):
                break
    return out


if __name__ == "__main__":
    print("d_s cross-check -- re-validating d_s-vs-d_H on CDT dual graphs\n")
    validate_against_dH()
    print("\nReminder: d_s is a QUALITATIVE cross-check here, not a verdict "
          "instrument. Quote the volume profile / cos^3 blob as the second "
          "observable; use d_s only to confirm it moves the same direction.")
