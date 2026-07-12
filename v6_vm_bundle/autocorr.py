#!/usr/bin/env python3
"""Autocorrelation-aware error estimation for the v6 analyzers.

Every quoted error in the paper's lattice sections currently carries the
caveat "not corrected for autocorrelation". This module removes it: the
integrated autocorrelation time tau_int is estimated with Sokal's adaptive
window on the FFT autocovariance, and standard errors are inflated by
sqrt(tau_int) (equivalently: n_eff = n / tau_int).

API:
  iat(x)                -> tau_int  (>= 1.0; 1.0 means uncorrelated)
  mean_err(x)           -> (mean, sem_corrected, tau_int, n_eff)
  combine_arms(means)   -> (mean, sem) across independent arms/seeds/pins
                           (plain scatter -- arms are independent by design)

Self-test: `python autocorr.py` checks the estimator on AR(1) chains with
known tau_int = (1+rho)/(1-rho) and on white noise.
"""
from __future__ import annotations
import numpy as np


def iat(x, c=6.0, max_win_frac=0.5):
    """Integrated autocorrelation time via Sokal's adaptive windowing.

    tau_int(W) = 1 + 2*sum_{t=1..W} rho(t), with the window W chosen
    self-consistently as the smallest W >= c*tau_int(W)  (Sokal, c~5-10).
    Returns 1.0 for series too short (or too flat) to estimate; never
    returns less than 1.0, so corrected errors are never smaller than the
    naive ones.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8:
        return 1.0
    x = x - x.mean()
    if not np.any(x):
        return 1.0
    # FFT autocovariance (zero-padded to avoid circular wrap)
    m = 1
    while m < 2 * n:
        m *= 2
    f = np.fft.rfft(x, m)
    acov = np.fft.irfft(f * np.conjugate(f))[:n].real
    if acov[0] <= 0:
        return 1.0
    rho = acov / acov[0]
    tau = 1.0
    wmax = max(2, int(n * max_win_frac))
    for W in range(1, wmax):
        tau = 1.0 + 2.0 * float(rho[1:W + 1].sum())
        tau = max(tau, 1.0)   # negative-rho noise can drive the sum < 1
        if W >= c * tau:
            break
    return float(max(1.0, tau))


def mean_err(x):
    """Mean with an autocorrelation-corrected standard error.

    Returns (mean, sem, tau_int, n_eff) where sem = std / sqrt(n_eff) and
    n_eff = n / tau_int. For n < 2 the sem is nan.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n == 0:
        return (float("nan"), float("nan"), float("nan"), 0.0)
    if n == 1:
        return (float(x[0]), float("nan"), 1.0, 1.0)
    tau = iat(x)
    neff = max(1.0, n / tau)
    sem = float(x.std(ddof=1) / np.sqrt(neff))
    return (float(x.mean()), sem, tau, float(neff))


def combine_arms(means):
    """Combine per-arm (per-seed / per-pin) means into a single value with a
    plain across-arm scatter error: mean +- std/sqrt(n_arms). Arms are
    independent by construction (different seeds or well-separated pins),
    so no autocorrelation correction applies at this level. This is the
    'quote across-seed scatter, not within-run SEM' rule of the runbook.
    Returns (mean, sem); sem is nan with fewer than 2 arms."""
    m = np.asarray([v for v in means if np.isfinite(v)], dtype=float)
    if m.size == 0:
        return (float("nan"), float("nan"))
    if m.size == 1:
        return (float(m[0]), float("nan"))
    return (float(m.mean()), float(m.std(ddof=1) / np.sqrt(m.size)))


def _selftest():
    rng = np.random.default_rng(0)
    print("autocorr self-test (AR(1): true tau_int = (1+rho)/(1-rho))")
    ok = True
    for rho, n in [(0.0, 20000), (0.5, 20000), (0.8, 40000), (0.9, 80000)]:
        true_tau = (1 + rho) / (1 - rho)
        x = np.empty(n)
        x[0] = rng.standard_normal()
        eta = rng.standard_normal(n)
        for i in range(1, n):
            x[i] = rho * x[i - 1] + eta[i]
        t = iat(x)
        rel = abs(t - true_tau) / true_tau
        good = rel < 0.15
        ok &= good
        print(f"  rho={rho:.1f} n={n}: tau_est={t:6.2f}  true={true_tau:6.2f} "
              f" rel.err={rel:.3f}  {'ok' if good else 'FAIL'}")
    # corrected sem must cover the true mean at the right rate
    rho, n, trials, hits = 0.8, 4000, 200, 0
    for _ in range(trials):
        x = np.empty(n)
        x[0] = rng.standard_normal()
        eta = rng.standard_normal(n)
        for i in range(1, n):
            x[i] = rho * x[i - 1] + eta[i]
        m, s, _, _ = mean_err(x)
        hits += abs(m) < 2 * s
    cov = hits / trials
    good = cov > 0.90
    ok &= good
    print(f"  2-sigma coverage at rho=0.8: {cov:.2f} (want ~0.95, > 0.90)"
          f"  {'ok' if good else 'FAIL'}")
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
