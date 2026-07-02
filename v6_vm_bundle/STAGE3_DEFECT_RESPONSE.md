# Stage 3 — defect response ("does geometry feel mass?")

The theory's definition of mass is persistent closure failure (committed
capacity; the electron formula reads mass as a closure-return rate of the
refresh dynamics). The coarse-graining results (`induced_couplings.py`,
`refresh_coarse_grain.py`) showed the closure sector cannot supply vacuum
curvature stiffness but *naturally* couples geometry to closure-failure
density. Stage 3 tests exactly that coupling: insert persistent failures
into a thermalized coupled universe and measure what happens around them.

## Design

Pin N well-separated spatial tetrahedra (`v6_defect_run.py`):

| arm | frozen labels | content |
|---|---|---|
| **real** (`--pin-mode fail`) | all four faces m=0 (6 collisions, K²=48) | the defect: maximal persistent closure failure |
| **placebo** (`--pin-mode closed`) | m={0,1,2,3} (injective, K²=40.67 = min) | identical anchoring + move-blocking, zero failure |
| **virtual** (measured inside both arms) | none | freshly sampled unpinned cells far from pins — the vacuum baseline |

Pinned triangles are excluded from the heat-bath; carrier tets are protected
two ways: (8,2) proposals on carrier vertices are vetoed pre-apply (vertex
ids are unstable under applied-then-undone (8,2), found the hard way), and
any other move that would remove a carrier is Metropolis-rejected via
dS=+1e30. **Pin survival is a hard gate: `pins alive: N/N` or the arm is
void.** Both arms carry the identical protection, so real-vs-placebo
subtracts the anchoring effect exactly.

Shell observables (BFS distance d = 1..rmax in the slice tet-adjacency),
accumulated every `--measure-every` sweeps:

- `mean_E` — mean closure energy of shell cells (**capacity strain**, the
  theory's gravity-analog field)
- `coll` — fraction of shell cells with a label collision
- `mean_q` — mean slice-edge coordination (**local curvature**)
- `shell_n` — cells per shell (**local volume**)

## Running it (per arm; run real + placebo as parallel processes)

```bash
# from a thermalized stage-2 beta=1 closure checkpoint:
python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode fail   \
       --pins 100 --target-n41 20000 --K 80 --sweeps 3000 --out def_real
python v6_defect_run.py --resume clo_b1.0_20k.json --pin-mode closed \
       --pins 100 --target-n41 20000 --K 80 --sweeps 3000 --out def_placebo
python v6_defect_run.py --analyze def_real.csv def_placebo.csv
```

~100 pins × ~1400 measurements each ≈ 10⁵ defect-environment samples per
arm; roughly overnight per arm at 20k (they can share the machine with
spare cores). The free preview (`failure_geometry_preview.py`, runs on
existing stage-2 checkpoints in seconds) gives a fluctuation-grade look at
the same question while you wait.

## Reading discipline — the claim ladder if it works

Grade the result by which observables separate (real vs placebo, decaying
toward virtual with distance):

1. **`mean_E` / `coll` separate** → *capacity response*: "persistent closure
   failure sources a strain field in the surrounding capacity state." This
   is the microscopic seed of the paper's source→strain mechanism —
   necessary for the gravity story, but it is a statement about the label
   sector, not yet about geometry.
2. **`mean_q` / `shell_n` also separate** → *geometric response*: "local
   geometry deforms specifically around the theory's mass, beyond the effect
   of anchoring alone." This is the headline-grade claim: the first
   dynamical demonstration that the theory's definition of matter couples to
   spacetime geometry.
3. **NOT available at any grade:** Newton's law, 1/r, G, or any long-range
   statement — the measured correlation length is sub-lattice, shells 1–2
   carry the signal, and these volumes cannot resolve a far field. Also not
   available: "mass curves spacetime" in the GR sense — the licensed sentence
   is "the defect deforms its local geometric environment."
4. **Null result** (nothing separates from placebo): also informative — the
   static closure coupling alone does not transmit failure into geometry at
   measurable strength, consistent with the ~1% induced-coupling estimate,
   and the geometry–matter story would then also need the non-local refresh
   structure (Many-Pasts), same as the vacuum sector.

Standard gates first, as always: foliation CLEAN, matched volume between
arms, pins alive N/N, audits passing. And the same phase caveat as stage 2:
until the phase grid lands, this runs on a near-boundary ensemble — a
detected response is meaningful (mechanism exists), a null is weak evidence.
