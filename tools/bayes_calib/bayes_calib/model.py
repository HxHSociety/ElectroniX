"""Forward model, prior, likelihood and posterior for Darveaux calibration.

The three ingredients of Bayes' rule, each in its own function:

    posterior(theta | data)  ∝  likelihood(data | theta)  ×  prior(theta)

theta = [ln K1, K2, ln K3, K4, ln beta]
        └────── Darveaux constants ──────┘  └ Weibull shape (scatter) ┘

K1 and K3 are positive and span orders of magnitude, so we infer their
LOGARITHMS. That keeps the sampler away from negative values and makes the
posterior closer to Gaussian, which the sampler handles best.

All functions are vectorised over a "walker" axis so the sampler can score
dozens of candidate parameter sets in one numpy call.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PARAM_NAMES = ["ln_K1", "K2", "ln_K3", "K4", "ln_beta"]
NATURAL_NAMES = ["K1", "K2", "K3", "K4", "beta"]

# Crack must grow across this fraction of the pad diameter to count as failure.
# MUST match rpim_solver/src/solder_fatigue.rs :: darveaux_life  (a_f = 0.65·d).
A_F_FRACTION = 0.65


# ─── Forward model ────────────────────────────────────────────────────────────


def ln_characteristic_life(theta: np.ndarray, dw: np.ndarray, d_mm: np.ndarray) -> np.ndarray:
    """ln(eta) = ln(Nf), the Weibull characteristic life, for every (walker, observation).

    Darveaux (same equations as the Rust solver):
        N0    = K1 · dW^K2                 crack-initiation cycles
        da/dN = K3 · dW^K4                 crack-growth rate [mm/cycle]
        Nf    = N0 + a_f / (da/dN)         with a_f = 0.65 · pad diameter

    Nf is a sum of two positive terms, so we work in logs with `logaddexp`
    to avoid overflow when the exponents are large.

    theta : (W, 5)   dw, d_mm : (N,)   ->   returns (W, N)
    """
    theta = np.atleast_2d(theta)
    ln_k1 = theta[:, [0]]
    k2 = theta[:, [1]]
    ln_k3 = theta[:, [2]]
    k4 = theta[:, [3]]

    ln_dw = np.log(np.maximum(dw, 1e-12))[None, :]
    ln_af = np.log(A_F_FRACTION * d_mm)[None, :]

    ln_n0 = ln_k1 + k2 * ln_dw                 # ln(K1 · dW^K2)
    ln_growth = ln_af - ln_k3 - k4 * ln_dw     # ln(a_f / (K3 · dW^K4))
    return np.logaddexp(ln_n0, ln_growth)


def ln_initiation_life(theta: np.ndarray, dw: np.ndarray) -> np.ndarray:
    """ln(N0): only the crack-initiation part (used to show identifiability)."""
    theta = np.atleast_2d(theta)
    return theta[:, [0]] + theta[:, [1]] * np.log(np.maximum(dw, 1e-12))[None, :]


# ─── Prior ────────────────────────────────────────────────────────────────────


@dataclass
class Prior:
    """Independent Gaussian priors + hard physical bounds.

    Defaults are CENTRED on the solver's built-in SAC305 constants
    (K1=22.18, K2=-1.26, K3=0.2169, K4=1.26, beta=5.4) with wide widths:
    "we believe the library is roughly right, but we are not sure by how much".
    The data then pull the posterior away from the prior only as far as the
    evidence demands.
    """

    mean: np.ndarray = None
    sd: np.ndarray = None

    def __post_init__(self) -> None:
        if self.mean is None:
            self.mean = np.array([np.log(22.18), -1.26, np.log(0.2169), 1.26, np.log(5.4)])
        if self.sd is None:
            self.sd = np.array([1.5, 0.5, 1.5, 0.5, 0.5])
        self.mean = np.asarray(self.mean, dtype=float)
        self.sd = np.asarray(self.sd, dtype=float)

    def log_prob(self, theta: np.ndarray) -> np.ndarray:
        theta = np.atleast_2d(theta)
        z = (theta - self.mean) / self.sd
        lp = -0.5 * np.sum(z * z, axis=1)
        # Physics: more strain energy => faster initiation (K2 < 0) and faster
        # crack growth (K4 > 0). Violations get zero prior probability.
        bad = (theta[:, 1] >= 0.0) | (theta[:, 3] <= 0.0)
        return np.where(bad, -np.inf, lp)


# ─── Likelihood ───────────────────────────────────────────────────────────────


@dataclass
class Dataset:
    """Fatigue-test observations. One row per specimen / daisy-chain.

    dw     : solver ΔW for that specimen's test condition [MPa]
    d_mm   : pad diameter [mm]
    cycles : cycles at failure OR at test suspension
    failed : True = failed at `cycles`; False = still alive (right-censored)
    group  : label of the test condition (for plots / reports)
    """

    dw: np.ndarray
    d_mm: np.ndarray
    cycles: np.ndarray
    failed: np.ndarray
    group: np.ndarray

    def __len__(self) -> int:
        return len(self.cycles)


def log_likelihood(theta: np.ndarray, data: Dataset) -> np.ndarray:
    """Weibull log-likelihood with right-censoring.

    A failed specimen contributes the Weibull density
        ln f(t) = ln β − β ln η + (β−1) ln t − (t/η)^β
    A suspended specimen (test stopped before failure) contributes only the
    probability of surviving that long, the survival function
        ln S(t) = −(t/η)^β
    Ignoring suspended units would bias the life estimate LOW, since only the
    early failures would be counted. Keeping them is the correct treatment.
    """
    theta = np.atleast_2d(theta)
    ln_beta = theta[:, [4]]
    beta = np.exp(ln_beta)

    ln_eta = ln_characteristic_life(theta, data.dw, data.d_mm)      # (W, N)
    ln_t = np.log(data.cycles)[None, :]                             # (1, N)

    z = np.exp(beta * (ln_t - ln_eta))                              # (t/η)^β
    ln_f = ln_beta - ln_eta + (beta - 1.0) * (ln_t - ln_eta)
    return np.sum(np.where(data.failed[None, :], ln_f, 0.0) - z, axis=1)


# ─── Posterior ────────────────────────────────────────────────────────────────


def make_log_posterior(data: Dataset, prior: Prior):
    """Return f(theta[W,5]) -> log posterior[W] (up to an additive constant)."""

    def log_posterior(theta: np.ndarray) -> np.ndarray:
        lp = prior.log_prob(theta)
        out = np.full(lp.shape, -np.inf)
        ok = np.isfinite(lp)
        if np.any(ok):
            ll = log_likelihood(theta[ok], data)
            out[ok] = lp[ok] + ll
        # NaN (e.g. overflow) must never be accepted by the sampler.
        return np.where(np.isnan(out), -np.inf, out)

    return log_posterior


def to_natural(theta: np.ndarray) -> np.ndarray:
    """[ln K1, K2, ln K3, K4, ln β]  ->  [K1, K2, K3, K4, β]."""
    theta = np.atleast_2d(theta).copy()
    out = theta.copy()
    out[:, 0] = np.exp(theta[:, 0])
    out[:, 2] = np.exp(theta[:, 2])
    out[:, 4] = np.exp(theta[:, 4])
    return out
