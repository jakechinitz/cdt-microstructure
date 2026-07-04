# The Regge bridge: a predicted volume-sector dressing (measured)

**Claim tested:** the substrate contributes a predicted extensive term to
the host Regge action's volume coupling. The closure sector's free energy
per cell (mu, computed analytically by transfer-matrix thermodynamic
integration from the paper's eta* and the injectivity penalty — no lattice
input) shifts the volume (cosmological/Delta) sector of the effective
Regge action. Any extensive matter free energy shifts a cosmological-sector
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
