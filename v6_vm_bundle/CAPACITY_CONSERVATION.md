# The lattice capacity-conservation law (theory → Regge, step 4)

**Question:** what conservation structure does the paper's finite-capacity
budget imply on the lattice — the fidelity anchor any future *propagation*
model (stage 4) must implement, playing the role g_share,eff = 7.4198 played
for the closure model.

**Answer in one line:** a budget alone gives *screened* (sub-cell-range)
response — no gravity; the paper's long-range sector requires a genuine
**continuity law** for free capacity plus its own **maintenance-flux**
reading of matter, and together those two produce the Newtonian form
Φ ∝ M/r structurally. Verified numerically on real slice geometry
(`capacity_lattice_law.py`).

## 1. What the paper fixes (inputs)

- **P1 — budget.** Finite capacity per cell; the vacuum cell's sharing
  entropy is the admissibility-closed g_share,eff = **7.4198** (re-derived
  by the script — the same anchor as the closure model).
- **P2 — commitment.** Matter locks capacity; committed capacity is
  conserved except at commitment/release events (the cosmology's
  commit/release phases).
- **P3 — maintenance.** A defect requires *continuous* maintenance, **linear
  in its commitment** ("the continuous maintenance of the defect is large
  and linear in the commitment while the rest energy … is exponentially
  small").
- **P4 — refresh bandwidth.** Every cell gets exactly one memoryless redraw
  per τ\* — an exactly conserved per-cell, per-tick resource.

## 2. Route A — budget bookkeeping alone: screened (no gravity)

Suppose free capacity is only *constrained*: link shares c_l ≥ 0 with
per-cell sums fixed, Σ_{l∋x} c_l = C − m_x, re-equilibrating by maximum
entropy. Stationarity gives s′(c_l) = μ_x + μ_y (one Lagrange multiplier per
cell); linearizing, δc_l = χ(δμ_x + δμ_y), and the cell constraint becomes

    χ [ z·δμ_x + Σ_{y~x} δμ_y ] = −m_x        (z = 4 faces/cell)
    ⇔  (2z·I − L) δμ = −m/χ                   (L = slice Laplacian)

The operator (2z − L) has **no small-k pole** — its Green function is
short-ranged with sub-cell decay. Numerically on a real slice (169 cells):
the response falls **eight orders of magnitude in 13 steps** (~×0.3 per
step). Conclusion: if the paper's capacity were merely a per-cell budget
that locally re-equilibrates, gravity would be exponentially screened at the
substrate scale. **Budget bookkeeping cannot be the paper's gravity.**

## 3. Route B — conservation + maintenance: graph-Coulomb (Newton's form)

Suppose instead free capacity obeys a genuine **continuity law** — it is
neither created nor destroyed, only transported through shared faces and
exchanged with commitment:

    f_x(t+dt) − f_x(t) = − Σ_{faces xy} J_xy − commit_x + release_x ,
    J_xy = −J_yx ,          Σ_x [ f_x + m_x ] = const  (exact)

and a defect draws a steady maintenance flux Q_x ∝ m_x (P3). The steady
state of transport (any local flux law linearizes to J ∝ −D∇f) solves

    L δf = (Q/D) (δ_source − uniform return)

— the **massless graph-Poisson equation**, massless *because* conservation
forbids a local restoring term. Its Green function on a 3D slice falls like
1/r (finite-size-flattened on our small closed slices; the numerical table
shows the contrast with Route A directly). The deficit field is then

    δf(r) ∝ Q/r ∝ m/r    ⇒    Φ ≡ δf/f̄ ∝ m/r ,

i.e. **Newton's form emerges structurally**, with the M-linearity of the
source inherited from P3 (maintenance linear in commitment) and the 1/r from
the conservation law — exactly the paper's weak-field bridge δS ↔ Φ.

## 4. The theorem-shaped conclusion

The paper's long-range gravity **requires the conservation reading**. This
is a selection result between two readings the paper's text supports, and it
tells us what the microscopic carrier must be. Two candidates:

- **B1 — capacity stock:** an amount of sharing capacity per cell,
  transported between neighbors. Needs a microscopic definition of the
  transported quantity beyond the label sector (the labels alone have no
  conserved quantity — measured fact from the refresh NESS analysis).
- **B2 — refresh bandwidth:** each cell's one-redraw-per-τ\* is *exactly*
  conserved per tick (P4); a defect's maintenance consumes nearby cells'
  refresh coherence, and the deficit field is the steady allocation of
  bandwidth. The paper leans this way: mass is a recurrence *rate*, a₀ ∝
  c·H(z) is rate-flavored, and maintenance-linear-in-commitment is natural
  for a bandwidth share. B2 needs no new conserved stuff — the conservation
  law is the tick itself.

Either way, the lattice law is the continuity equation of §3 with the
following **fidelity gates for any stage-4 propagation model**:

- **T1 (exactness):** Σ_x [f_x + m_x] constant to machine precision under
  the dynamics, changing only at commit/release events.
- **T2 (vacuum anchor):** homogeneous steady state with ⟨f⟩ = g_share,eff =
  7.4198 per cell.
- **T3 (Coulomb check):** a pinned steady sink of strength Q produces a
  deficit profile matching the graph-Poisson Green function (1/r on large
  slices), with amplitude linear in Q. This is the long-range analog of
  stage 3's near-field measurement, and stage 3's pinned-defect machinery is
  the natural host for the sink.
- **T4 (screening diagnostic):** any leakage — capacity created or destroyed
  locally — shows up as Yukawa screening of the T3 profile. Fit a screening
  mass; it must be consistent with zero. This is the model's version of the
  placebo discipline: it detects the failure mode (Route A contamination)
  rather than assuming its absence.

## 5. Status and honest scope

This is a *derivation of requirements*, not yet a model: it identifies the
conserved structure the paper needs, shows the alternative dies by
screening, and hands stage 4 its gates. It does not (and cannot) decide
B1 vs B2 — that is a theory-side choice, though the paper's own language
favors B2. And it changes nothing about the earlier no-gos: the vacuum's
geometric stiffness remains an input; what conservation buys is the
*long-range reach of the matter-sourced strain field*, which is the sector
the paper actually claims. Chain so far: closure owns what mass is; refresh
owns how mass strains its neighborhood (stage 3); conservation owns how the
strain travels (stage 4); vacuum geometry is an input at every step.
