# Running the v6 4D CDT test on a VM

This bundle tests whether the v6 Causal-Dynamical-Triangulations engine generates
**4-dimensional spacetime** from the graph (the "substrate" / emergence claim),
first with bare Einstein-Regge dynamics (the known-answer gate) and then with the
theory's EPRL spin-foam amplitude on top.

> **Why this is a VOLUME test (and not a strictness test).** We checked locally:
> at strictness levels 0, 1, and 2 the engine gives *identical* geometry and stays
> perfectly simplicial — strictness is **not** the lever. The one hypothesis left
> for "why only ~2D so far" is **volume**: the de Sitter (4D) phase of CDT is a
> large-N phenomenon (literature onset N_4 ≳ 60k, production runs 160k-360k, with
> ~80 time slices). So the experiment below is a **volume scan**: grow N_41 and
> watch whether the Hausdorff dimension climbs toward the 4-torus rail and a
> localized de Sitter "blob" forms.

---

## 0. Prerequisites

* Python **3.10+** (tested on 3.12).
* `pip install numpy scipy`  (nothing else — the engine is pure Python/NumPy).
* A machine with **3+ cores** is ideal (run the volume scan in parallel).

## 1. The bundle (13 files — all required)

```
v6_cdt.py  v6_cdt_moves.py  v6_cdt_run.py  v6_run_lib.py
v6_verify_run.py  v6_theory_run.py
v5_s1xs3_init.py  v5_cdt_link.py  v5_cdt_lib.py  v5_cdt_moves.py
step3_linkA_harness.py  vertex_tensor.py  vertex_j3.npz
```

Unzip into one directory and `cd` into it. Quick smoke test (≈30 s):

```bash
python -u v6_verify_run.py --target-n41 200 --K 8 --max-sweeps 3 --measure-every 3
# expect: "manifold check ... PASS" and a VERDICT line.
```

---

## 2. THE MAIN TEST — bare-Regge volume scan (no theory yet)

This is the gate. Run three volumes **in parallel** (one process per core). Each
checkpoints to disk and prints a live table; `nohup ... &` lets them run for
hours/days and survive disconnects.

```bash
# 10k, 20k, 40k target N_41, all at K=80 time slices, canonical de Sitter point.
nohup python -u v6_verify_run.py --target-n41 10000 --K 80 --k0 2.2 --Delta 0.6 \
      --k4 0.9 --eps 1e-4 --checkpoint scan_10k.json  > scan_10k.log 2>&1 &
nohup python -u v6_verify_run.py --target-n41 20000 --K 80 --k0 2.2 --Delta 0.6 \
      --k4 0.9 --eps 1e-4 --checkpoint scan_20k.json  > scan_20k.log 2>&1 &
nohup python -u v6_verify_run.py --target-n41 40000 --K 80 --k0 2.2 --Delta 0.6 \
      --k4 0.9 --eps 1e-4 --checkpoint scan_40k.json  > scan_40k.log 2>&1 &
```

Watch progress:  `tail -f scan_20k.log`

### Live table columns

```
 sweep    N4    N41   d_H   blob  active  cos3err  valid  links   wall_s
```

* **N41** should climb to the target and then hover (the `eps` penalty pins it).
* **d_H** — Hausdorff dimension of the dual graph. THE headline number.
* **blob** — peak/mean of the spatial-volume profile (1 = flat stalk, >>1 = blob).
* **active** — # time slices clearly above the stalk (de Sitter extent).
* **cos3err** — relative error of an A·cos³ fit to the profile (lower = more de-Sitter).
* **valid / links** — must read `ok` / `ok` (or `gen`); a `BAD` means a bug — stop.

### How to read the final VERDICT (printed at the end / on resume-and-finish)

The estimator under-reads, so judge **relative to the rails** it prints:

```
rails (matched size): 2-torus≈1.9   3-torus≈2.6   4-torus≈3.1
```

* **Success (graph makes 4D):** as target N_41 increases 10k→20k→40k, **d_H climbs
  monotonically toward the 4-torus rail**, the **blob score rises**, **active
  slices grow**, and **cos3err falls**. That trend *is* the finite-size-scaling
  signature of emergent de Sitter space — the result you want.
* **Still collapsed:** d_H stays near the 2-torus rail and the profile stays flat
  (`active` ≈ 1-2, blob ≈ 1) at all three volumes. Then either even 40k is too
  small (push higher / longer) or the engine genuinely doesn't reach phase C.

Always confirm the final line:
`manifold check (gluing-based + S^3 links): PASS`  — if it says FAIL, the geometry
isn't a valid manifold and the d_H number is meaningless.

### Resuming a run (after a stop / disconnect / to grow further)

```bash
python -u v6_verify_run.py --resume scan_20k.json --checkpoint scan_20k.json \
       --max-sweeps 1000000 > scan_20k.cont.log 2>&1 &
```

You can also cap wall time:  add `--wall-hours 12`.

---

## 3. THE THEORY TEST — only meaningful AFTER step 2 shows 4D

If (and only if) the bare gate reaches 4D, run the EPRL amplitude **on top** and
compare to the bare run **at the same volume**. The robust signal is the
*difference* (finite-size effects cancel): does the amplitude steer the geometry?

> **Post-audit protocol (2026-07).** The coupling machinery was audited and
> fixed — see `SIM_AUDIT_coupling_misspecifications.md`. The sampler is now
> exact for a declared joint measure: β-consistent label heat-bath,
> free-energy centering (mu via thermodynamic integration, printed as
> `# [mu TI]` lines and stored in the checkpoint), slot-symmetrized vertex
> tensor, and a **mandatory shuffled-tensor placebo arm**. Use the sweep
> script — it wires all of this up:

```bash
# from the bundle directory, with the bare checkpoint from step 2 available:
BASE=scan_20k.json ./run_eprl_sweep.sh "0.0 0.05 0.1 0.3" 20000
# watch:
tail -n 4 logs/thy_b*_20k.log
grep -hE 'mu TI|centering|audit ok|Traceback|BAD' logs/thy_b*_20k.log | tail -n 30
```

* Each arm resumes from the SAME thermalized bare checkpoint, pins the same
  N41, calibrates its own mu (a few minutes of heat-bath passes at startup),
  and runs with `--local-eprl --audit-every 25` (recompute-and-abort-on-drift
  safety net).
* The placebo arm (`thy_b0.3shuf_20k`) runs the largest β with an
  entry-shuffled tensor (`make_shuffled_control.py`): same entry statistics,
  all EPRL structure destroyed.

### How to read it (gates BEFORE observables)

1. **G3′ — matched N4.** All arms must sit at comparable **N4** (not just the
   pinned N41). If N4 trends with β, the term is not volume-neutral — stop
   and investigate; do not read d_H.
2. **Placebo separation.** Compare the real β_max arm to the shuffled arm.
   Only a real-vs-shuffled **difference** is attributable to the EPRL
   amplitude's structure; if they match, the movement (or the flatness) is
   machinery/entry-statistics, not EPRL.
3. Then read d_H / blob / cos3err vs β against the bare run, as before.

### Honest caveats (printed by the theory run, repeated here)

1. The vertex amplitude is the **frozen-j=3** tensor (`vertex_j3.npz`),
   positivized (|A| — the true intertwiner contraction has order-unity sign
   interference) and slot-symmetrized (the raw tensor's slot order is an
   unchecked convention). It is a PLACEHOLDER for the EPRL amplitude;
   **absolute** numbers are artefact-grade. What's meaningful is the
   matched-volume comparison read against the placebo.
2. A frozen-j amplitude has no area degrees of freedom (all triangles carry
   the same spin), so it cannot feel deficit angles; the measured intertwiner
   variance is type-blind. A flat sweep is therefore expected for structural
   reasons and does NOT vindicate or falsify the theory (see the audit doc's
   frozen-j vs peaked-j section).
3. `--k4 0.9` is a starting guess in our conventions. If N_41 runs away or
   freezes, retune `--k4` (and/or `--eps`); the volume penalty is the real
   anchor.

---

## 4. Rough sizing

Pure-Python, single-threaded. The prior v6 run reached N_4≈6400 in ~30 min.
Expect roughly: 10k in ~1-2 h, 20k in a few h, 40k in many h. The literature's
60k onset may be slow here — but the **trend across 10k/20k/40k** is what answers
the question, even short of full 4D. Run the three volumes on separate cores.

To self-check the engine before a long run (≈2 min): `python -u v6_cdt_moves.py`
— runs all 7 Pachner-move inverse/manifold tests + an all-moves stress test
(should end with `all moves preserve a valid CDT 4-manifold = True`).
