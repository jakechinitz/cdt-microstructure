# ============================================================================
#  STEP 3, LINK A: VALIDATED SPECTRAL-DIMENSION HARNESS
#
#  This is the MEASUREMENT INSTRUMENT for the manifoldlikeness test. Its whole
#  purpose is to be TRUSTWORTHY: before it measures anything novel, it proves it
#  can recover KNOWN spectral dimensions and -- the decisive capability -- cleanly
#  distinguish a branched polymer (d_s = 4/3, the theory-killer) from a
#  4-manifold (d_s = 4, the theory-confirmer).
#
#  WHY THIS EXISTS: a novel d_s number has no answer key. So we manufacture one:
#  we calibrate the estimator's (systematic, monotonic) finite-size bias against
#  tori of KNOWN dimension, then VALIDATE that the calibration generalizes to a
#  held-out structure (branched polymer) it was not fit on. Only a harness that
#  passes this earns the right to measure the unknown condensate.
#
#  VALIDATED RESULT (this file, as run):
#    calibration  d_s_true ~ 0.877 * d_s_raw - 0.007   (fit on tori d=2,3,4)
#    held-out branched polymer: calibrated 1.28 vs known 1.33  -> generalizes
#    discrimination: killer 1.28 vs confirmer 3.98, separation 2.70 -> CLEAN
#
#  Link B (the condensate model that feeds this harness) is built separately and
#  measured with THIS calibration. Do not measure a novel model with an
#  uncalibrated estimator.
# ============================================================================

import numpy as np
from scipy.sparse import csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply

# ---------------------------------------------------------------- estimator ---
def raw_ds(adj, tmin=0.5, tmax=40, nt=30, n_probe=80, seed=0):
    """Raw d_s(t) by sparse heat-kernel return probability (stochastic trace).
       Scales to large sparse graphs (no dense eigendecomposition)."""
    rng = np.random.default_rng(seed)
    A = csr_matrix(adj.astype(float)); N = A.shape[0]
    deg = np.asarray(A.sum(1)).ravel(); deg[deg==0] = 1; s = 1/np.sqrt(deg)
    L = diags(np.ones(N)) - diags(s) @ A @ diags(s)
    ts = np.logspace(np.log10(tmin), np.log10(tmax), nt)
    pr = rng.standard_normal((N, n_probe)) / np.sqrt(n_probe)
    P = np.array([np.mean(np.sum(pr * expm_multiply(-t*L, pr), 0)) for t in ts])
    return ts, -2*np.gradient(np.log(np.clip(P,1e-12,None)), np.log(ts))

def midwin(ts, ds, lo=3, hi=15):
    w = (ts>lo) & (ts<hi); return float(np.median(ds[w]))

# ------------------------------------------------------ known test structures --
def torus(side, dim):
    N = side**dim; grid = np.arange(N).reshape((side,)*dim); A = lil_matrix((N,N))
    for ax in range(dim):
        r = np.roll(grid, 1, axis=ax)
        for i,j in zip(grid.ravel(), r.ravel()): A[i,j]=A[j,i]=1
    return A.tocsr()

def branched_polymer(N, seed=0):
    rng = np.random.default_rng(seed); A = lil_matrix((N,N))
    for i in range(1,N):
        j = rng.integers(0,i); A[i,j]=A[j,i]=1
    for _ in range(N//20):
        a,b = rng.integers(0,N,2)
        if a!=b: A[a,b]=A[b,a]=1
    return A.tocsr()

# ----------------------------------------------------------- calibration -------
def fit_calibration(verbose=True):
    """Fit d_s_true ~ a*d_s_raw + b on tori of known dimension. Returns (a,b)."""
    known = np.array([2.0, 3.0, 4.0])
    raw = []
    for side,dim in [(24,2),(11,3),(6,4)]:
        ts,ds = raw_ds(torus(side,dim)); raw.append(midwin(ts,ds))
    raw = np.array(raw)
    M = np.vstack([raw, np.ones_like(raw)]).T
    a,b = np.linalg.lstsq(M, known, rcond=None)[0]
    if verbose: print(f"calibration: d_s_true ~ {a:.3f}*d_s_raw + {b:.3f}")
    return a, b

def validate(a, b):
    """Confirm the calibration recovers tori AND the held-out branched polymer,
       and that killer vs confirmer are cleanly separated."""
    cal = lambda r: a*r + b
    print(f"\n{'structure':26}{'calibrated':>12}{'known':>8}{'verdict':>9}")
    allpass = True
    checks = [("2-torus",torus(24,2),2.0),("3-torus",torus(11,3),3.0),
              ("4-torus",torus(6,4),4.0),("branched polymer*",branched_polymer(2000),4/3)]
    bp = m4 = None
    for name,G,known in checks:
        ts,ds = raw_ds(G); c = cal(midwin(ts,ds))
        ok = abs(c-known) < 0.35; allpass &= ok
        if "branched" in name: bp = c
        if "4-torus" in name: m4 = c
        print(f"{name:26}{c:>12.2f}{known:>8.2f}{'PASS' if ok else 'CHECK':>9}")
    print(f"\nDISCRIMINATION: killer {bp:.2f} vs confirmer {m4:.2f} "
          f"-> separation {m4-bp:.2f} ({'CLEAN' if m4-bp>2 else 'TOO CLOSE'})")
    print(f"HARNESS: {'VALIDATED -- trustworthy for novel measurement' if allpass else 'NEEDS WORK'}")
    return allpass

def measure(adj, a, b, label="model"):
    """The trusted measurement: calibrated d_s of a novel structure."""
    ts, ds = raw_ds(adj)
    c = a*midwin(ts,ds) + b
    print(f"  [{label}] calibrated d_s = {c:.2f}  "
          f"-> {'4-MANIFOLD (theory lives)' if c>3.3 else ('BRANCHED POLYMER (theory dies)' if c<2 else 'INTERMEDIATE/unclear')}")
    return ts, ds, c

if __name__ == "__main__":
    print("STEP 3 LINK A -- validating the d_s measurement instrument\n")
    a, b = fit_calibration()
    ok = validate(a, b)
    print("\nLink A provides: raw_ds(), fit_calibration(), measure(adj,a,b).")
    print("Link B (the condensate model) is measured by calling measure() with")
    print("THIS (a,b). A number from measure() is only trustworthy because this")
    print("harness proved the instrument recovers known answers first.")
