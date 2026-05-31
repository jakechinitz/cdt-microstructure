#!/usr/bin/env python3
"""EPRL term diagnostic: is the (centered) vertex-amplitude term INERT or does
it have TEETH?

The matched-volume beta sweep can only be interpreted if the centered EPRL term
actually carries variance across configurations. If the centered
per-pentachoron log-amplitude is a near-delta spike, the term is inert.

>>> CORRECTION (do not over-read "HAS TEETH") <<<
Spread alone (sigma>0) is NOT sufficient to make a flat sweep interpretable. What
matters is whether the variance couples to the GEOMETRIC degrees of freedom (the
(4,1)/(3,2) simplex-type structure) that d_H and the volume profile are built
from. Measured on vertex_j3.npz:
    sigma (per-pent centered)         = 0.87   (real spread)
    type-split |<41> - <32>|          = 0.054  (16x SMALLER than sigma)
    dual-graph nearest-nbr corr r_nn  = +0.14  (weak, ~4 sigma vs shuffled null)
So the variance is LARGELY (not perfectly) orthogonal to geometry. Consequence:
a FLAT d_H sweep is WEAK evidence -- flat is the expected outcome almost
regardless of beta when the coupling rides mostly along label directions d_H
cannot see. Do NOT bank a flat sweep as "theory tolerates 4D"; that is hard to
distinguish from "the term cannot move d_H." The informative question is whether
the weak-but-real geometric channel moves the cos^3 BLOB (not just d_H) as beta
rises to 0.3.

Verified: this type-blindness is NOT an artifact of centering -- subtracting a
global scalar mu cannot change a difference of means; raw and centered
type-splits are identical to 4e-15. The orthogonality is a property of the
frozen-j=3 amplitude itself.

This is a ZERO-SWEEP test: it builds one thermalized-label config and histograms
  centered_p = -log|A_v(faces of p)| - mu ,   mu = mean over pentachora.

Reported:
  * sigma = std of centered contributions
  * the effective per-move action scale beta*sigma at the swept betas
  * the (4,1)/(3,2) type-split -- if << sigma, the teeth are largely orthogonal
    to the geometry d_H measures, so a flat sweep is weak (not strong) evidence.

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
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--betas", default="0.05 0.1 0.3")
    args = ap.parse_args()

    Tt = FaithfulVertex.load(args.vertex).dense_tensor().astype(np.float64)
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
    hb = make_heatbath(intw, Tt)
    for _ in range(args.hb):
        hb(T, rng)

    raw, typ = per_pent_neglog(T, intw, Tt)
    mu = raw.mean()
    cen = raw - mu
    sigma = cen.std()
    print(f"\n[heat-bath equilibrated, N4={len(raw)}]  mu={mu:.4f}")
    print(f"  centered -log|amp|: sigma={sigma:.4f}  min={cen.min():.3f} max={cen.max():.3f}")
    qs = np.percentile(cen, [1, 10, 25, 50, 75, 90, 99])
    print(f"  percentiles [1,10,25,50,75,90,99] = {np.round(qs, 3)}")

    intw_u = Intertwiners(T, D, rng)
    raw_u, _ = per_pent_neglog(T, intw_u, Tt)
    print(f"  (uniform-random labels sigma = {(raw_u-raw_u.mean()).std():.4f}  "
          f"-> heat-bath narrows the distribution)")

    # raw (uncentered) type-split too: proves centering does not create blindness
    raw41, raw32 = raw[typ == "41"], raw[typ == "32"]
    m41, m32 = cen[typ == "41"], cen[typ == "32"]
    if len(m41) and len(m32):
        split_cen = abs(m41.mean() - m32.mean())
        split_raw = abs(raw41.mean() - raw32.mean())
        print(f"\n[geometry coupling] mean centered by simplex type: "
              f"(4,1)={m41.mean():+.4f}  (3,2)={m32.mean():+.4f}")
        print(f"  type-split: RAW={split_raw:.4f}  CENTERED={split_cen:.4f}  "
              f"(diff={abs(split_raw-split_cen):.1e} -> centering does NOT cause blindness)")
        print(f"  split/sigma = {split_cen/sigma:.3f}  "
              f"-> {'type-BLIND: teeth largely orthogonal to geometry; a flat sweep is WEAK evidence' if split_cen < 0.2*sigma else 'type-COUPLED: term distinguishes the geometric DOF'}")

    # spatial structure on the dual graph (the geometry d_H actually traverses)
    pids = list(T.pent.keys()); idx = {p: i for i, p in enumerate(pids)}
    xs, ys = [], []
    for p in pids:
        for q in T.nbr[p]:
            if q in idx and idx[q] > idx[p]:
                xs.append(cen[idx[p]]); ys.append(cen[idx[q]])
    if xs:
        r_nn = float(np.corrcoef(xs, ys)[0, 1])
        rng2 = np.random.default_rng(1)
        null = [np.corrcoef(cen, rng2.permutation(cen))[0, 1] for _ in range(200)]
        sig_n = float(np.std(null))
        print(f"\n[dual-graph structure] nearest-nbr corr r_nn={r_nn:+.4f}  "
              f"(null std {sig_n:.4f})  -> "
              f"{'WHITE NOISE (geometry-blind)' if abs(r_nn) < 3*sig_n else f'weak structure ({abs(r_nn)/sig_n:.0f} sigma): real but small geometric channel'}")

    print(f"\n[verdict] sigma={sigma:.4f} (spread present), but geometric coupling "
          f"is WEAK (split/sigma small, r_nn small).")
    print("  => A flat d_H sweep is WEAK evidence. Watch the cos^3 BLOB, not just")
    print("     d_H, for whether the small geometric channel bites as beta->0.3.")
    for b in (float(x) for x in args.betas.split()):
        print(f"  beta={b}: per-pent action scale beta*sigma = {b*sigma:.4f}")


if __name__ == "__main__":
    main()
