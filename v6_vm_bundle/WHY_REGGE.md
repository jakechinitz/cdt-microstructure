# Why Regge — the architecture of the theory, as computed

This is the synthesis document for the theory→lattice program: what the
simulation campaign and the six analytic calculations established about the
relationship between the entropic-substrate theory and the background
(Regge) geometry it is tested on. It is the reasoning a reader needs before
interpreting ANY run in this repo, and the skeleton of the paper's
computational-companion section.

## 1. Two different things called "spacetime emergence"

**Kinematic emergence — the substrate ENCODES geometry, by construction.**
Every cell is a geometrically consistent quantum tetrahedron: the
admissibility closure K² is precisely the closure condition of the quantum
tetrahedron, and the paper's length anchor λₑ/L\* = (2/3)e^{7g} converts
entanglement increments into meters. IMPORTANT NUANCE (paper-faithful, and
load-bearing for Lorentz invariance): the cells are DIMENSIONLESS
combinatorial objects — units of capacity, not space-voxels with rigid
edges. "The granularity lives in the capacity, not in a preferred spatial
lattice" (the paper's own phrasing), which is exactly what evades the
photon/GW dispersion bounds that killed rigid-lattice discrete-spacetime
models. L\* is a statistical conversion rate, not a grid pitch. The
simulation respects this separation: the theory's cells are the label
configurations; the lattice simplices merely HOST them (supplying
adjacency), as regulator scaffolding — and all measurements are matched
comparisons at fixed regulator, so scaffold-scale artifacts cancel. This
encoded-geometry claim cannot fail and no result in this repo questions it.

**Dynamical selection — why is the realized geometry a smooth, extended,
4D de Sitter one** rather than any of the combinatorially overwhelming
crumpled/branched alternatives? This is a *ranking* job. A ruler cannot
rank. Something must make the fabric lie flat.

Confusing these two produces both false alarm ("the theory can't make
spacetime!") and false comfort ("spacetime is built in!"). Every result
below is about the second job only.

## 2. The six-branch no-go: nothing specified does the ranking

Each branch is a computed result with a script in this repo:

| # | candidate ranking mechanism | result | where |
|---|---|---|---|
| 1 | static closure weight (η\*, injectivity) | induced curvature coupling ~1% of need; correlation range < 1 cell | `induced_couplings.py` |
| 2 | Markovian refresh (Model R NESS) | vacuum steady state geometry-blind to exponential accuracy | `refresh_coarse_grain.py` |
| 3 | closure-weighted histories (Many-Pasts tilt, Doob-exact) | link correlation capped ~0.13; structural ceiling at any strength | `many_pasts_requirements.py` |
| 4 | free conserved capacity field, back-coupled | induced action = ½·log det′L; per-cell value near-universal at fixed degree (Δ ≈ 2·10⁻⁴) | `field_induced_action.py` |
| 5 | budget-bounded (non-Gaussian) field | saturation mass m² = 2/g *suppresses* the only geometry-sensitive modes; closes tighter | `nongaussian_field_check.py` |
| 6 | (necessity theorem) local label dynamics cannot even equilibrate one cell | 48 frozen sectors of 35, verified exactly | `many_pasts_fracture_check.py` |

Uncovered containers, named: long-range field kernels and non-equilibrium
geometry–field co-evolution — both unspecified by the paper; testing an
invented version would test the invention (the EPRL-proxy lesson).

## 3. The paper's own answer — and why it makes Regge legitimate

The paper never claims its dynamics select the vacuum. §26.8, explicitly:
the framework *"does not sum over arbitrary triangulations,"* and *"the
history weighting is supplied by Many-Pasts rather than by a path integral
over geometries."* The vacuum's 4D form enters by **conditioning on the
present** — we observably inhabit an extended 4D universe, and the theory's
coefficients are computed within that conditioned state.

**Regge calculus is the standard mathematical spelling of "we are in a
smooth 4D universe."** Running the substrate on a Regge/CDT background is
therefore not the adoption of foreign physics — it is the only Monte Carlo
implementation of the paper's own Postulate-III conditioning. And the
six-branch no-go upgrades this from authorial choice to computed necessity:
the division of labor (conditioning supplies the vacuum; the substrate
supplies matter and response) is *forced* by the theory's specified
ingredients. The paper's architecture and the mathematics agree,
independently.

## 4. The relation is not decorative: the substrate dresses Regge at a predicted size

Three rungs, honestly graded:

1. **MEASURED — the volume-sector dressing** (`REGGE_BRIDGE.md`). The
   closure sector's free energy per cell, computed analytically from η\*
   alone, contributes an extensive term to the host action's volume
   coupling, with magnitude predicted in advance: ΔN41 = −β·μ/(4ε) = −62.2.
   (Any extensive matter free energy shifts a cosmological-sector coupling;
   the content is the predicted magnitude, not the shift's existence.) Three-arm measurement:
   uncentered − centered = **−62.8** (~1%). A UV-spec number appears as a
   gravitational-action coefficient in a dynamical measurement.
2. **PREDICTED — curvature-sector shifts.** c₀ ≈ 0.019 (closure) and
   ~0.08/cell (7-channel field, local part) are computed; the measurement
   is a phase-boundary displacement between bare and coupled ensembles at
   production statistics. Pre-registered, not yet resolved.
3. **STRUCTURAL — kinship, not identity.** K² measures the failure of a
   cell's oriented faces to close; a Regge deficit angle measures the
   failure of the dihedral angles around a hinge to close. The kinship is
   structural rather than identical — an honest analogy, and the smallness
   of rung 2 is the measured statement of how weakly the two levels
   currently couple.

## 5. What the substrate demonstrably does, inside the conditioned vacuum

- **Mass is closure failure** — the model reproduces the paper's numbers
  exactly (⟨K²⟩ = 3/(2η\*), g = 7.4198) as fidelity gates.
- **Mass strains its neighborhood** — stage 3, placebo-controlled, both
  capacity and local-geometry observables (single-seed caveats recorded).
- **A conserved medium carries the strain to distance** — stage 4;
  conservation exact; commitment ladder monotone; source law and junction
  constant one v1.1 run from claim-grade; Coulomb/screening gates (T3/T4)
  owed at volume. The Coulomb/Yukawa dichotomy (`CAPACITY_CONSERVATION.md`)
  shows Newton's *form* Φ ∝ m/r assembles from conservation + the
  maintenance postulate.
- **The medium does nothing to empty space** — provably, branches 1–5.
  Gravity in this framework is response around matter, never vacuum
  self-organization. This matches the paper's stated scope.

## 6. The claim discipline (what any write-up may say)

Demonstrated: mechanism results with placebo separation, inside the
conditioned vacuum. Computed: every no-go branch and the bridge
prediction. Gated/open: the engine's own de Sitter calibration (phase
grid, autotuned, in progress — a tooling matter, not a referendum on the
theory, per §3); rung-2 coupling shifts; stage-4 T3/T4; replication seeds.
Never available from this program: the value of G, galaxy-scale
predictions, or the quantum-interpretation content of Postulate III.