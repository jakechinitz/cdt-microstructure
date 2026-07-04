# The Regge bridge: a predicted volume-sector dressing (measured)

**Claim tested:** the substrate contributes a predicted extensive term to
the host Regge action's volume coupling. Integrating the labels out at
fixed geometry defines S_eff[g] = S_host[g] - log Z_label[g]; the extensive
part of -log Z_label is beta*mu(beta)*N_cells, with mu the per-cell label
free energy (computed by annealed thermodynamic integration on the base
geometry before any coupled dynamics run; the exact zero-lattice single-cell
value from direct 7^4 enumeration, 2.477 from eta* and the injectivity
penalty alone, is within 0.5%). That term shifts the volume
(cosmological/Delta) sector of the effective Regge action. Any extensive matter free energy shifts a cosmological-sector
coupling — the content here is not the existence of the shift but its
magnitude, predicted in advance:
coupling the substrate uncentered at beta=1 must displace the equilibrium
volume, relative to the centered arm, by

    delta(N41) = -beta * mu / (4 * eps)
               = -2.4883 / 0.04  =  -62.2      (demo settings: eps = 0.01)

**Measurement** (three arms, identical protocol, common clean base at
N41 = 2000, K = 16, 150 sweeps, late-run means):

| arm                  | mean N41 | shortfall |
|----------------------|----------|-----------|
| bare (beta = 0)      | 1993.2   |  -6.8     |
| centered (beta = 1)  | 1969.8   | -30.2     |
| uncentered (beta = 1)| 1907.0   | -93.0     |

    uncentered - centered = -62.8    vs    predicted -62.2   (~1%)

The subtraction isolates exactly the mu-term (both arms share all other
machinery). The centered arm's own -30 is the known second-order
non-neutrality of the TI centering (real, gated by matched-N4 in
production; it cancels in the difference).

**Reading:** a number computed from the theory's UV spec (eta*, the
seven-state ensemble) appears, at the predicted magnitude, as a coefficient
of the Regge action's volume sector in a dynamical lattice measurement.
The substrate does not derive Regge — the six-branch no-go stands — but it
demonstrably *dresses* Regge: a computable, predicted-in-advance piece of
the volume coupling (this measurement), and computable local shifts of the curvature-
sector couplings (c0 ~ 0.019 closure, ~0.08 seven-channel field; phase-
boundary-displacement measurement pre-registered for production
statistics). Regge-as-background is therefore not an arbitrary import:
it is the conditioned vacuum of the paper's Many-Pasts architecture, and
the substrate detectably renormalizes it in the sectors the theory owns.

Reproduce: three v6_closure_run.py arms from one clean checkpoint --
beta 0 / beta 1 centered / beta 1 --no-center-closure -- and compare
late-run mean N41; prediction from calibrate_mu_ti at the same settings.

## The scaling family (pre-registered)

The demonstration point is one member of a predicted family,

    Delta N41(beta, eps) = -beta * mu(beta) / (4 * eps),

with three separately testable signatures:

- **S1, beta-curve (concavity + Ehrenfest).** beta*mu(beta) =
  int_0^beta <E>_s ds and the Gibbs mean falls as labels order, so
  Delta(beta) is concave — e.g. at eps=0.01 the prediction is -37 at
  beta=0.5, not half of the -61 at beta=1. Local slope:
  d(Delta)/d(beta) = -<E>_beta/(4*eps), where <E>_beta is measurable in
  the same arm whose displacement it predicts.
- **S2, eps-law.** Delta is exactly proportional to 1/eps at fixed beta.
  Departures measure the curvature of the background free energy that the
  matched-pair subtraction cancels at first order — the derivation's one
  assumption, now itself instrumented.
- **S3, placebo tracks its own table.** The orbit-shuffled placebo table
  has its own free energy: single-cell mu_plc = 3.17 vs real mu = 2.48
  (28% apart). The placebo pair is predicted to displace by
  -beta*mu_plc/(4*eps) — not by zero, and not by the real value. The
  number must track the table, which upgrades the placebo from a null
  control to a second quantitative point.

Tooling: `bridge_predict.py` (one TI pass per table gives the whole
mu(beta) curve; writes the prediction CSV before any arm launches) and
`run_bridge_sweep.sh` (matched centered/uncentered pairs per (beta, eps)
point from one clean base, placebo pair at the largest beta; analyze with
`bridge_predict.py --analyze logs/bridge_*.log --pred <csv>`). Gates:
prediction/measurement ratio ~ 1 across the family, concave beta-curve,
1/eps collapse, placebo on the mu_plc line, foliation CLEAN on every arm.
