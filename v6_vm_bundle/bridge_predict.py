#!/usr/bin/env python3
"""Regge-bridge scaling family: predictions (default) and log analysis.

THE FAMILY. Integrating the closure labels out at fixed geometry defines an
effective action for the host, S_eff[g] = S_host[g] - log Z_label[g]. The
extensive part of -log Z_label is beta*mu(beta)*N_cells (mu = per-cell label
free energy); the centered ensemble subtracts it, the uncentered ensemble
keeps it, and against the volume pin eps*(N41 - Nbar)^2 the difference must
displace the equilibrium volume by

    Delta N41(beta, eps) = -beta*mu(beta) / (4*eps).

The demonstration point (beta=1, eps=0.01: predicted -62.2, measured -62.8)
is one member. Three pre-registered signatures upgrade it to a family test:

  S1  beta-curve: Delta(beta) is CONCAVE, since beta*mu(beta) =
      int_0^beta <E>_s ds and the Gibbs mean falls as labels order; local
      slope d(Delta)/d(beta) = -<E>_beta/(4*eps), an Ehrenfest relation in
      which the same run supplies both sides (<E>_beta is measurable in the
      arm whose displacement it predicts).
  S2  eps-law: Delta is exactly proportional to 1/eps at fixed beta;
      departures measure the curvature of the background free energy that
      the matched-pair subtraction cancels at first order (the derivation's
      one assumption, now itself instrumented).
  S3  placebo tracks its own table: the orbit-shuffled placebo table has
      its own free energy mu_plc (single-cell ~3.17 vs real ~2.48 at the
      theory point, 28% apart), so the placebo arm is predicted to displace
      by -beta*mu_plc/(4*eps) -- not by zero and not by the real value.
      The measured number must track the table it was computed from.

PREDICT mode computes mu(beta) for the real and placebo tables by ONE
annealed thermodynamic-integration pass each on the base geometry
(cumulative trapezoid gives the whole beta-curve), prints the predicted
displacement grid, and writes a CSV for the analyzer. The exact zero-lattice
single-cell value (direct 7^4 enumeration; 2.477 at the theory point, within
~0.5% of the lattice value) is printed as an analytic reference.

ANALYZE mode (--analyze LOG...) parses logs from run_bridge_sweep.sh
(tags bridge_b{B}_e{E}_{cen|unc}[plc]_*), averages N41 over the last --tail
fraction of measurement rows, pairs centered/uncentered arms, and tables
measured vs predicted displacement using the CSV from a predict pass.

Usage:
  python bridge_predict.py --resume scan_20000_causal.json \
      --betas "0.25 0.5 1.0" --eps-list "5e-4 1e-3 2e-3" --out bridge_pred.csv
  python bridge_predict.py --analyze logs/bridge_*.log --pred bridge_pred.csv
"""
from __future__ import annotations
import argparse
import re
import numpy as np

from v6_closure_run import (build_energy_table, ETA_STAR, FaceLabels,
                            tet_of, tris_of, _heatbath_pass)

_trapz = getattr(np, "trapezoid", None) or np.trapz


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def single_cell_mu(Etab, betas, n=400):
    """Exact zero-lattice reference: TI over the 7^4 single-cell ensemble."""
    E = Etab.reshape(-1).astype(float)
    out = {}
    for b in betas:
        b = float(b)
        if b <= 0:
            out[b] = (float(E.mean()), float(E.mean()))
            continue
        ss = np.linspace(0.0, b, n)
        means = []
        for s in ss:
            w = np.exp(-s * (E - E.min()))
            means.append(float((E * w).sum() / w.sum()))
        out[b] = (float(_trapz(means, ss) / b), means[-1])
    return out


def ti_curve(T, Etab, betas, seed, subdiv=4, equil=3, measure=2, verbose=True):
    """One annealed TI pass over s in [0, max(betas)] on the geometry T.
    Returns {beta: (mu(beta), <E>_beta)} at every requested beta (cumulative
    trapezoid: beta*mu(beta) = int_0^beta <E>_s ds), and the tet count."""
    rng = np.random.default_rng(seed)
    labels = FaceLabels(rng)
    tet_pids = {}
    for p, vs in T.pent.items():
        tet = tet_of(T, vs)
        if tet is not None:
            tet_pids.setdefault(tet, set()).add(p)
    for tet in tet_pids:
        labels.ensure(tris_of(tet))

    def mean_E():
        s = 0.0
        for tet in tet_pids:
            l = [labels.lab[tri] for tri in tris_of(tet)]
            s += Etab[l[0], l[1], l[2], l[3]]
        return s / max(1, len(tet_pids))

    bs = sorted(set(float(b) for b in betas if float(b) > 0))
    grid = [0.0]
    prev = 0.0
    for b in bs:
        grid.extend(np.linspace(prev, b, int(subdiv) + 1)[1:].tolist())
        prev = b
    grid = np.array(grid)
    fs = []
    for s in grid:
        for _ in range(int(equil)):
            _heatbath_pass(labels, tet_pids, Etab, s, rng)
        acc = []
        for _ in range(int(measure)):
            _heatbath_pass(labels, tet_pids, Etab, s, rng)
            acc.append(mean_E())
        fs.append(sum(acc) / len(acc))
        if verbose:
            print(f"#   s={s:.4f}  <E>_s = {fs[-1]:.4f}", flush=True)
    fs = np.array(fs)
    cum = np.concatenate([[0.0],
                          np.cumsum((fs[1:] + fs[:-1]) / 2 * np.diff(grid))])
    out = {}
    for b in betas:
        b = float(b)
        if b <= 0:
            out[b] = (float(fs[0]), float(fs[0]))
            continue
        i = int(np.argmin(np.abs(grid - b)))
        out[b] = (float(cum[i] / grid[i]), float(fs[i]))
    return out, len(tet_pids)


def predict(args):
    betas = [float(b) for b in args.betas.split()]
    eps_list = [float(e) for e in args.eps_list.split()]
    if args.resume:
        from v6_run_lib import load_checkpoint
        T, _, _ = load_checkpoint(args.resume)
        geom = f"{args.resume} (N4={T.n_pent()})"
    else:
        from v6_cdt import build_s1xs3
        T = build_s1xs3(K=args.K)
        geom = f"fresh thin S1xS3, K={args.K} (production: use the sweep's base checkpoint)"
    print(f"# Regge-bridge predictor  eta={args.eta:.6g}  lambda_inj={args.lambda_inj}")
    print(f"# base geometry: {geom}")

    tables = [("real", build_energy_table(args.eta, args.lambda_inj)),
              ("placebo", build_energy_table(args.eta, args.lambda_inj,
                                             args.placebo_seed))]
    rows = []
    for name, Etab in tables:
        print(f"# [{name}] annealed TI pass (subdiv={args.subdiv}, "
              f"equil={args.equil}, measure={args.measure}) ...")
        latt, ntet = ti_curve(T, Etab, betas, seed=args.seed + 777,
                              subdiv=args.subdiv, equil=args.equil,
                              measure=args.measure, verbose=args.verbose)
        exact = single_cell_mu(Etab, betas)
        print(f"# [{name}] {ntet} spatial tets; zero-lattice single-cell "
              f"reference in parentheses")
        print(f"#   {'beta':>6} {'mu_latt':>9} {'(mu_1cell)':>11} {'<E>_beta':>9}")
        for b in betas:
            mu, Eb = latt[float(b)]
            mu1, _ = exact[float(b)]
            print(f"#   {b:>6.3g} {mu:>9.4f} {mu1:>11.4f} {Eb:>9.4f}")
        for b in betas:
            mu, Eb = latt[float(b)]
            for e in eps_list:
                rows.append((name, float(b), float(e), mu, Eb,
                             -float(b) * mu / (4 * e),
                             -Eb / (4 * e)))
    print(f"\n# predicted uncentered-minus-centered displacement "
          f"Delta N41 = -beta*mu(beta)/(4*eps)")
    print(f"{'table':>8} {'beta':>6} {'eps':>9} {'mu':>9} {'DeltaN41':>10} "
          f"{'dDelta/dbeta':>13}")
    for r in rows:
        print(f"{r[0]:>8} {r[1]:>6.3g} {r[2]:>9.3g} {r[3]:>9.4f} "
              f"{r[5]:>10.2f} {r[6]:>13.2f}")
    if args.out:
        with open(args.out, "w") as f:
            f.write("table,beta,eps,mu,E_beta,delta_pred,slope_pred\n")
            for r in rows:
                f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.6f},{r[4]:.6f},"
                        f"{r[5]:.4f},{r[6]:.4f}\n")
        print(f"# wrote {args.out}")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

_TAG = re.compile(r"bridge_b([0-9.]+)_e([0-9.eE+-]+)_(cen|unc)(plc)?")
_ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s")
_MU_TI = re.compile(r"mu\([0-9.eE+-]+\)\s*=\s*(-?[0-9.]+)")
_MU_REUSE = re.compile(r"reusing mu\s*=\s*(-?[0-9.]+)")


def logged_mu(path):
    """The mu the centered arm actually subtracted (last TI/reuse line in its
    log). The pair difference isolates exactly this constant -- whatever the
    TI's own noise -- so -beta*mu_arm/(4*eps) is the EXACT prediction for the
    pair, while the predictor CSV is the independent pre-registration."""
    mu = None
    with open(path) as f:
        for ln in f:
            m = _MU_TI.search(ln) or _MU_REUSE.search(ln)
            if m:
                mu = float(m.group(1))
    return mu


def tail_mean_n41(path, tail):
    vals = []
    with open(path) as f:
        for ln in f:
            m = _ROW.match(ln)
            if m:
                vals.append(int(m.group(3)))
    if not vals:
        return None, 0
    k = max(1, int(len(vals) * tail))
    return float(np.mean(vals[-k:])), k


def analyze(args):
    pred = {}
    if args.pred:
        with open(args.pred) as f:
            next(f)
            for ln in f:
                t, b, e, mu, Eb, d, sl = ln.strip().split(",")
                pred[(t, float(b), float(e))] = float(d)
    arms = {}
    for path in args.analyze:
        m = _TAG.search(path)
        if not m:
            print(f"# skip (no bridge tag): {path}")
            continue
        b, e = float(m.group(1)), float(m.group(2))
        arm, plc = m.group(3), bool(m.group(4))
        mean, k = tail_mean_n41(path, args.tail)
        if mean is None:
            print(f"# skip (no measurement rows): {path}")
            continue
        arms.setdefault((b, e, plc), {})[arm] = (mean, k, path)
    print(f"{'table':>8} {'beta':>6} {'eps':>9} {'N41_cen':>9} {'N41_unc':>9} "
          f"{'measured':>9} {'pre-reg':>9} {'exact':>9} {'ratio':>6}")
    for (b, e, plc), d in sorted(arms.items()):
        if "cen" not in d or "unc" not in d:
            print(f"# incomplete pair beta={b} eps={e} plc={plc}: "
                  f"have {sorted(d)}")
            continue
        meas = d["unc"][0] - d["cen"][0]
        table = "placebo" if plc else "real"
        p = pred.get((table, b, e))
        mu_arm = logged_mu(d["cen"][2])
        px = (-b * mu_arm / (4 * e)) if mu_arm is not None else None
        ref = px if px is not None else p
        ratio = (meas / ref) if ref else float("nan")
        print(f"{table:>8} {b:>6.3g} {e:>9.3g} {d['cen'][0]:>9.1f} "
              f"{d['unc'][0]:>9.1f} {meas:>9.1f} "
              f"{(f'{p:.1f}' if p is not None else '--'):>9} "
              f"{(f'{px:.1f}' if px is not None else '--'):>9} {ratio:>6.2f}")
    print("# 'pre-reg' = predictor-CSV value (independent TI, on record "
          "before launch); 'exact' = -beta*mu_arm/(4*eps) using the mu the "
          "centered arm actually subtracted (its own log) -- the pair "
          "difference isolates exactly that constant, so ratio uses it.")
    print("# gates: ratio ~ 1 across the family; concave beta-curve (S1); "
          "1/eps collapse (S2); placebo tracks mu_plc, not mu_real (S3).")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--betas", type=str, default="0.25 0.5 1.0")
    p.add_argument("--eps-list", type=str, default="5e-4 1e-3 2e-3")
    p.add_argument("--eta", type=float, default=ETA_STAR)
    p.add_argument("--lambda-inj", type=float, default=3.0)
    p.add_argument("--placebo-seed", type=int, default=12345)
    p.add_argument("--resume", type=str, default=None,
                   help="base checkpoint (the one the sweep arms resume from)")
    p.add_argument("--K", type=int, default=16)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--subdiv", type=int, default=6,
                   help="TI grid points per beta interval")
    p.add_argument("--equil", type=int, default=3)
    p.add_argument("--measure", type=int, default=2)
    p.add_argument("--out", type=str, default="bridge_pred.csv")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--analyze", nargs="*", default=None,
                   help="run-log files to analyze instead of predicting")
    p.add_argument("--pred", type=str, default="bridge_pred.csv",
                   help="prediction CSV from a predict pass (analyze mode)")
    p.add_argument("--tail", type=float, default=0.5,
                   help="fraction of late measurement rows to average")
    args = p.parse_args()
    if args.analyze:
        analyze(args)
    else:
        predict(args)


if __name__ == "__main__":
    main()
