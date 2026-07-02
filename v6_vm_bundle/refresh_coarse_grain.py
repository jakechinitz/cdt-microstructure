#!/usr/bin/env python3
"""FROM THEORY TOWARD REGGE, step 2: coarse-grain the REFRESH dynamics.

Step 1 (induced_couplings.py) showed the STATIC closure weight induces a
curvature-sector coupling ~1% of what phase C needs. The theory's remaining
dynamical ingredient is Postulate III's memoryless refresh: each cell, at
Poisson rate 1/tau*, redraws its whole boundary as a fresh independent draw
from the admissibility ensemble p_eta*(b) ~ exp(-eta* K^2) over the 1680
injective states (the paper's App. D.4 shows local moves cannot do this --
the refresh is irreducibly non-local at the cell scale).

MODEL R (the declared fork): on a slice complex, cells refresh independently;
a SHARED face carries the value written by whichever adjacent cell refreshed
most recently ('last writer wins'). This is the unconditioned reading of
'each pass is an independent draw'. The alternative reading -- redraw
CONDITIONED on the current neighbor faces -- is a heat-bath block sampler
whose stationary state is exactly the static Gibbs measure, i.e. it
coarse-grains back to step 1 and adds nothing new. Model R is therefore the
only reading that could add dynamics beyond the static weight; it produces a
genuine nonequilibrium steady state (NESS).

EXACT NESS STRUCTURE (single cell). Let the cell's age since its own last
refresh be a ~ Exp(1) (rates in refresh units). Each face is 'kept' iff its
neighbor has not refreshed within that age, independently with prob e^{-a}.
Integrating the age out, the probability that a SPECIFIC set of k faces is
kept (the rest overwritten) is exactly

    P_k = k! (4-k)! / 5!        (sums correctly: sum_k C(4,k) P_k = 1)

Kept faces carry the cell's own correlated p_eta draw; overwritten faces
carry i.i.d. copies of the single-face marginal q1 of p_eta (the neighbor's
draw is exchangeable, and its refresh time tells us nothing about its
labels). Everything below follows by enumeration over the 840 injective
tuples -- no simulation, no approximation beyond Model R itself.

COMPUTED OBSERVABLES:
  * P_inj^NESS  -- probability a vacuum cell is currently ADMISSIBLE
                   (injective). Fresh cells are admissible by construction;
                   the NESS value is < 1: the vacuum carries a finite
                   density of TRANSIENT closure failures. This is the
                   model's vacuum defect density -- the 'mass-like' order
                   parameter of empty space, and the natural quantitative
                   home for 'mass is the probability of closure(-return)':
                   persistent (pinned) failure = matter; transient failure =
                   vacuum fluctuation.
  * <K^2>^NESS  -- mean closure defect in the NESS vs the ensemble value
                   3/(2 eta*) = 50.223.

GEOMETRY-BLINDNESS (the structural result). In Model R the vacuum NESS is
independent of the slice's gluing pattern to exponential accuracy:
  1. every cell has exactly 4 face-neighbors regardless of geometry, so the
     keep/overwrite law P_k is the same everywhere;
  2. overwriting values are i.i.d. marginals -- they carry no information
     about WHICH neighbor wrote them, hence none about the graph;
  3. coordination (q) sensitivity requires correlations transmitted AROUND
     a ring; transmission needs the same neighbor to have last-written two
     ring faces AND its draw to correlate them -- each link costs a factor
     (writer overlap ~ 1/3) x (within-draw pair correlation), so ring terms
     scale like c^q with c << 1, far below even step 1's induced coupling.
Therefore: vacuum refresh in Model R induces NO curvature-sector action.
What it DOES induce, for any geometry dynamics gated on local closure
(moves easier where cells close), is a coupling of geometry to the local
CLOSURE-FAILURE DENSITY -- i.e. to matter. The refresh dynamics naturally
generates the gravity-RESPONSE sector (geometry feels mass), not the vacuum
Regge stiffness. This matches the paper's own architecture: its gravity is
capacity strain AROUND committed defects (a response), and it explicitly
lacks -- and does not claim -- a microscopic vacuum action (S26.8 open
task). The two no-go branches together (conditioned refresh -> static
weight -> ~1%; unconditioned refresh -> geometry-blind vacuum) say the
vacuum 4D stiffness must be an input at this level of the theory, with the
refresh sector owning matter and geometry-matter coupling.

Usage:  python refresh_coarse_grain.py [--eta ETA_STAR]
"""
from __future__ import annotations
import argparse
from itertools import permutations, product
from math import factorial
import numpy as np

from v6_closure_run import ETA_STAR


def k2(m):
    S = sum(m)
    Sig2 = sum(x * x for x in m)
    return 48.0 - (S * S - Sig2) / 3.0


def ness_report(eta=ETA_STAR, verbose=True):
    # the admissibility ensemble over ordered injective 4-tuples
    tuples = list(permutations(range(-3, 4), 4))          # 840
    w = np.array([np.exp(-eta * k2(t)) for t in tuples])
    w /= w.sum()
    K2s = np.array([k2(t) for t in tuples])

    # single-face marginal q1 (exchangeable -> position irrelevant)
    q1 = np.zeros(7)
    for t, p in zip(tuples, w):
        q1[t[0] + 3] += p

    # exact keep-set law: P(specific k faces kept) = k!(4-k)!/120
    Pk = {k: factorial(k) * factorial(4 - k) / 120.0 for k in range(5)}

    P_inj = 0.0
    K2_mean = 0.0
    for t, p in zip(tuples, w):
        for keep in product((0, 1), repeat=4):
            k = sum(keep)
            kept = [t[i] for i in range(4) if keep[i]]
            n_over = 4 - k
            pk = Pk[k]
            # enumerate overwrites (i.i.d. q1) over the free coordinates
            for xs in product(range(-3, 4), repeat=n_over):
                px = 1.0
                for x in xs:
                    px *= q1[x + 3]
                cfg = kept + list(xs)
                weight = p * pk * px
                if len(set(cfg)) == 4:
                    P_inj += weight
                # reconstruct full config in original positions for K^2
                # (K^2 is symmetric, order irrelevant)
                K2_mean += weight * k2(cfg)

    if verbose:
        print(f"[Model R NESS at eta={eta:.6g}]")
        print(f"  single-face marginal q1 = {np.round(q1, 4)}")
        print(f"  P(cell admissible) in NESS      = {P_inj:.4f}")
        print(f"  vacuum transient-failure density = {1 - P_inj:.4f}")
        print(f"  <K^2> in NESS = {K2_mean:.3f}   vs ensemble 3/(2 eta) = "
              f"{3 / (2 * eta):.3f}")
        print(f"  (fresh cells: admissible by construction; the NESS deficit"
              f" is the refresh dynamics' intrinsic vacuum defect density)")
    return {"P_inj": P_inj, "fail": 1 - P_inj, "K2": K2_mean, "q1": q1}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eta", type=float, default=ETA_STAR)
    args = ap.parse_args()
    ness_report(args.eta)
