# Vertex amplitude provenance — `vertex_j3.npz`

This document records how the EPRL vertex amplitude tensor shipped in this repo
(`v6_vm_bundle/vertex_j3.npz`) was actually produced, what was and was not
validated, and the known reproducibility gaps. It exists because the amplitude
is the physical heart of the theory test, and the generating code does **not**
live in this repository — so without this record the tensor is a black box.

> **Read this before quoting any EPRL result.** The honest caveats below
> (especially §3 and §5) bound what the theory runs can legitimately claim.

---

## 1. What the file is

`vertex_j3.npz` is a NumPy archive holding the rank-5 EPRL vertex amplitude for
a 4-simplex, as a dense tensor over the five tetrahedral intertwiner labels.

| Field | Value | Meaning |
|---|---|---|
| `tensor` | shape `(7,7,7,7,7)`, float64 | `A_v(i0,i1,i2,i3,i4)` over 5 intertwiners |
| `j` | `3` | **frozen** boundary spin on every triangle |
| `immirzi` | `1.2` | Barbero–Immirzi parameter γ (field-standard) |
| `shells` | `5` | booster truncation Dl used in production |
| `y_map` | `rho_gj` | EPRL Y-map convention (γj boost) |
| `validated` | `1` | the small-spin reliability gate passed (see §3) |

Intertwiner dimension D = 7 = 2j+1. Tensor `max|A_v|` ≈ 8.77e-12, inside the
documented ~1e-13–1e-12 plausibility band for j=3 at γ=1.2.

The amplitude is computed with **sl2cfoam-next**
(`https://github.com/qg-cpt-marseille/sl2cfoam-next`).

---

## 2. How it was actually generated (the real path)

There are **two** code paths in the original project (`C:\Paper Sim\_archive\sl2cfoam_pipeline\`).
The production tensor came through the second one. This distinction matters:

- **Narrative / documented path:** `sl2cfoam_driver.jl` — a clean single-file
  Julia driver. **This is NOT what produced the shipped tensor.** Its hardcoded
  `SHELLS = 3` and its `j^12` asymptotic gate do **not** match the artifact
  (which is `shells=5` and was gated differently — see §3).

- **Actual production path:**
  1. sl2cfoam-next emits a `.sl2t` binary (the raw amplitude).
  2. `sl2cfoam_to_npz.py` converts `.sl2t` → `vertex_j3.npz`.
  3. `sl2cfoam_validate.py` runs the reliability gate and stamps `validated`.

When reproducing or extending this, follow the `.sl2t → to_npz → validate`
route, not the Julia driver's narrative.

Generation parameters (from the stored npz metadata, authoritative over any
hardcoded values in the driver): **j = 3, γ = 1.2, Dl (shells) = 5,
Y-map = rho_gj**.

---

## 3. What "validated = 1" means — and what it does NOT mean

This is the single most important caveat in this document.

`validated = 1` means the **small-spin reliability gate passed.** That gate, in
`sl2cfoam_validate.py`, is a 3-part check:
1. **Booster convergence** — Dl=3 vs Dl=5 stability of the amplitude.
2. **Structural sanity** — correct shape, finite, real.
3. **Magnitude plausibility** — `max|A_v|` ~ 1e-13–1e-12 for j=3 at γ=1.2.

`validated = 1` does **NOT** mean the asymptotic Regge gate passed. The classic
large-spin check (`|A_v| · j^12 → const`, EPRL → Regge asymptotics) was
**deliberately not run**, and correctly so: per `sl2cfoam_validate.py`, that gate
is invalid at small spin —

> *"That gate is wrong for small spin: EPRL Regge asymptotics kick in only at
> large j (typically j > 20). Trying it at j=3..5 fails for physics reasons, not
> library bugs."*

The large-spin asymptotic correctness of sl2cfoam-next is therefore **trusted
from the maintainers' published validation, not reproduced in-house.** That is a
defensible choice, but for any external claim it is a *stated assumption*, not an
in-house verification.

**Convergence evidence preserved on disk** (the booster-convergence test, part 1):
| File | shells (Dl) | validated | role |
|---|---|---|---|
| `vertex_j3_dl1.npz` | 1 | 0 | under-converged, correctly **rejected** |
| `vertex_j3.npz` (shipped) | 5 | 1 | converged, **accepted** |

The Dl=1-vs-Dl=5 pair *is* the reproducible convergence signal. Sibling tensors
`vertex_j2.npz` (Dl=5) and `vertex_j4.npz` (Dl=3) also exist.

`vertex_validation.log` (the driver's intended human-readable log) was **not
preserved.** Only the binary `validated` flag survived; the underlying
flatness/convergence numbers are gone.

---

## 4. Known reproducibility gap: unpinned library version

**No sl2cfoam-next commit, tag, or version string is recorded anywhere.** The
only reference is the package URL in the driver header. Exact bit-for-bit
reproduction therefore requires reconstructing the original install environment,
which is not captured in this repo. The tensor itself — a `(7,7,7,7,7)` float64
array plus the five metadata fields in §1 — is the only durable artifact.

If exact reproducibility becomes necessary, the missing pieces to recover are:
the sl2cfoam-next commit hash, the `.sl2t` source binary, and the
`sl2cfoam_to_npz.py` / `sl2cfoam_validate.py` scripts (currently in
`C:\Paper Sim\_archive\sl2cfoam_pipeline\`, not in this repo).

---

## 5. Frozen-j vs peaked-j: what this tensor can and cannot test

`vertex_j3.npz` is the **frozen-j = 3 placeholder**, not a peaked-j amplitude.
Every boundary triangle carries the same spin j = 3; only the five intertwiner
labels vary. The full theory's behavior is expected to depend on a *peaked-j*
superposition over multiple spins — and **those multi-j tensors were never
built.** Only single-j tensors (j = 2, 3, 4) exist.

Empirical behavior of this exact tensor (from in-session shell tests of the v5
engine, recorded for pre-registration):
1. It couples **strongly** to geometry (raw per-move ΔS_EPRL ~ 14–20× ΔS_Regge;
   note this is the *uncentered* scale — the v6 runs use a mean-subtracted /
   centered action, which greatly reduces the effective coupling).
2. At frozen j = 3 it tends to steer geometry toward **branched / sparse**
   configurations — i.e. d_H / d_s *downward*, away from 4.
3. This downward push is the **expected frozen-j behavior**, NOT a falsification
   of the theory. The theory's own escape hatch — that a peaked-j amplitude
   corrects the direction — requires the multi-j tensors that do not yet exist.

### Pre-registered interpretation (lock this in before reading results)

At matched volume, with a thermalized and unfrozen Markov chain:

| Observation | Interpretation |
|---|---|
| EPRL pushes d_H/d_s **down** (branched/sparse) | **Expected frozen-j behavior.** Consistent with theory + known peaked-j gap. **Does NOT falsify.** |
| d_H/d_s stays **flat** at the bare value | EPRL fluctuations don't steer at frozen-j. Also consistent. |
| EPRL pushes d_H/d_s **up** toward 4 | A *surprise* — frozen-j is not expected to help. |
| Markov chain frozen (acceptance ~0) | **Not a result** — coupling too strong; discard. |

**Consequence for scope:** a frozen-j run can only *check that frozen-j behaves
as the peaked-j argument predicts*. It **cannot falsify the theory**, because the
frozen-j truncation is a known, expected-to-be-wrong-direction placeholder. A
genuine falsification test requires building the peaked-j (multi-j) amplitude
first. Treat any frozen-j result as a consistency check, not a verdict.

---

## 6. Companion validation architecture (Step 3 pipeline)

A separate, more self-disciplined test architecture exists (see
`PIPELINE_README.md` if committed): the "Step 3 / Link B condensate" pipeline,
which measures the *spectral* dimension d_s through a validated estimator
(`step3_linkA_harness.py`) and gates its verdict behind three pre-registered
trust gates (faithful vertex / seed variance > 0.1 / thermalization). That
pipeline shares this same `vertex_j3.npz` but is a different model from the v6
CDT engine currently run on the VM (which measures the Hausdorff dimension d_H).
The two are complementary; neither inherits the other's trust gates
automatically.

---

*Maintained alongside the code so the honest caveats travel with the tensor
instead of living in a chat log. If you regenerate or extend the amplitude,
update §1–§4; if you build peaked-j tensors, revisit §5.*
