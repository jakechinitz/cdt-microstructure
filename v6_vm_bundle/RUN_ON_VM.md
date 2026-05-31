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

```bash
# matched-volume comparison at, say, N_41 = 20000:
#  (a) bare reference is scan_20k.json from step 2.
#  (b) Regge + a small EPRL correction (safer first test), fast local action:
nohup python -u v6_theory_run.py --mode regge_plus_eprl --beta-eprl 0.1 \
      --local-eprl --audit-every 25 --target-n41 20000 --K 80 \
      --checkpoint thy_b0.1_20k.json > thy_b0.1_20k.log 2>&1 &
```

* `--mode regge_plus_eprl` adds `beta_eprl * S_EPRL` to the Regge action (safe).
  Sweep `--beta-eprl` 0.0 → 0.1 → 0.3 → 1.0 and watch d_H / blob / cos3err move.
* `--local-eprl` uses the O(footprint) incremental action; `--audit-every 25`
  recomputes from scratch every 25 sweeps and **aborts on any drift** (safety net).
  Drop `--local-eprl` to use the simple global O(N) action as a cross-check.
* `--mode eprl_only` replaces Regge entirely (only interpretable once the amplitude
  fidelity is validated — see the caveats printed at startup).

### Honest caveats (printed by the theory run, repeated here)

1. The vertex amplitude is the **frozen-j=3** tensor (`vertex_j3.npz`); the full
   {15j} contraction still needs convention-checking, so **absolute** d_H from the
   theory run is artefact-grade. What's robust is the **bare-vs-EPRL comparison at
   matched volume**.
2. New tetrahedra get a fresh (uniform) intertwiner not folded into the Hastings
   ratio; the per-sweep heat-bath re-equilibrates labels. Fine for the
   steer-direction comparison, not for precision sampling.
3. `--k4 0.9` is a starting guess in our conventions. If N_41 runs away or freezes,
   retune `--k4` (and/or `--eps`); the volume penalty is the real anchor.

---

## 4. Rough sizing

Pure-Python, single-threaded. The prior v6 run reached N_4≈6400 in ~30 min.
Expect roughly: 10k in ~1-2 h, 20k in a few h, 40k in many h. The literature's
60k onset may be slow here — but the **trend across 10k/20k/40k** is what answers
the question, even short of full 4D. Run the three volumes on separate cores.

To self-check the engine before a long run (≈2 min): `python -u v6_cdt_moves.py`
— runs all 7 Pachner-move inverse/manifold tests + an all-moves stress test
(should end with `all moves preserve a valid CDT 4-manifold = True`).
