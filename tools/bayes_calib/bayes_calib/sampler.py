"""Affine-invariant ensemble MCMC sampler (Goodman & Weare, 2010).

Why MCMC at all?
----------------
The posterior p(theta | data) has no closed form: the Darveaux life is a sum
of two power laws inside a Weibull likelihood. MCMC draws *samples* from it
instead. Once you hold, say, 20 000 samples of theta, every question becomes
an average: the mean of K1 is the mean of the K1 column, a 90 % interval is
two percentiles, and the uncertainty of a prediction is that prediction
evaluated at each sample.

Why this particular sampler?
----------------------------
K1/K3 and K2/K4 are strongly correlated, so the posterior is a long thin
ridge. A plain random-walk Metropolis sampler crawls along it. The ensemble
"stretch move" is invariant to linear rescaling: it proposes new points for
one walker along the line joining it to another walker, so the proposal
shape adapts to the ridge automatically, with no tuning matrix.

The algorithm (one sweep)
-------------------------
Split walkers into halves A and B. To update every walker in A (in parallel):
  1. pick a random partner Xj from B
  2. draw a stretch factor z with density g(z) ∝ 1/sqrt(z) on [1/a, a]
  3. propose  Y = Xj + z (Xk − Xj)
  4. accept with probability min(1, z^(D−1) · p(Y) / p(Xk))
Then swap the roles of A and B. The z^(D−1) factor is what keeps the move
reversible in D dimensions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class Chain:
    samples: np.ndarray       # (steps, walkers, dim)
    log_prob: np.ndarray      # (steps, walkers)
    acceptance: float         # overall acceptance fraction


def ensemble_sample(
    log_prob: Callable[[np.ndarray], np.ndarray],
    p0: np.ndarray,
    n_steps: int,
    a: float = 2.0,
    seed: int | None = None,
) -> Chain:
    """Run the stretch-move ensemble sampler.

    log_prob : maps (W, D) -> (W,) log posterior density (vectorised)
    p0       : (W, D) starting positions; W must be even and >= 2*D
    """
    rng = np.random.default_rng(seed)
    n_walkers, dim = p0.shape
    if n_walkers % 2 or n_walkers < 2 * dim:
        raise ValueError("need an even number of walkers, at least 2 x dimension")

    cur = p0.astype(float).copy()
    cur_lp = log_prob(cur)
    if not np.all(np.isfinite(cur_lp)):
        raise ValueError("initial positions must all have finite log probability")

    half = n_walkers // 2
    halves = (slice(0, half), slice(half, n_walkers))

    samples = np.empty((n_steps, n_walkers, dim))
    lps = np.empty((n_steps, n_walkers))
    n_accept = 0

    for step in range(n_steps):
        for active, other in (halves, halves[::-1]):
            partners = cur[other][rng.integers(0, half, half)]
            # z ~ g(z) ∝ 1/sqrt(z) on [1/a, a], by inverse-CDF sampling
            z = ((a - 1.0) * rng.random(half) + 1.0) ** 2 / a
            proposal = partners + z[:, None] * (cur[active] - partners)
            prop_lp = log_prob(proposal)

            log_ratio = (dim - 1) * np.log(z) + prop_lp - cur_lp[active]
            accept = np.log(rng.random(half)) < log_ratio   # NaN compares False

            block = cur[active]
            block_lp = cur_lp[active]
            block[accept] = proposal[accept]
            block_lp[accept] = prop_lp[accept]
            cur[active] = block
            cur_lp[active] = block_lp
            n_accept += int(accept.sum())

        samples[step] = cur
        lps[step] = cur_lp

    return Chain(samples, lps, n_accept / (n_steps * n_walkers))


# ─── Diagnostics ──────────────────────────────────────────────────────────────


def _acf(x: np.ndarray) -> np.ndarray:
    """Normalised autocorrelation function of a 1-D series via FFT."""
    n = len(x)
    x = x - x.mean()
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    return acf / acf[0] if acf[0] > 0 else np.zeros(n)


def autocorr_time(samples: np.ndarray, c: float = 5.0) -> np.ndarray:
    """Integrated autocorrelation time per parameter [steps].

    MCMC samples are correlated: consecutive draws are not independent. The
    autocorrelation time tau says how many steps make ONE independent draw,
    so the effective sample size is  ESS = (steps x walkers) / tau.
    Uses Sokal's automatic windowing (stop summing once M >= c·tau).
    """
    n_steps, _, dim = samples.shape
    taus = np.empty(dim)
    for j in range(dim):
        acf = np.mean([_acf(samples[:, w, j]) for w in range(samples.shape[1])], axis=0)
        tau_running = 2.0 * np.cumsum(acf) - 1.0
        window = np.arange(len(tau_running)) >= c * tau_running
        m = int(np.argmax(window)) if window.any() else len(tau_running) - 1
        taus[j] = tau_running[m]
    return taus


def effective_sample_size(samples: np.ndarray, tau: np.ndarray) -> np.ndarray:
    n_steps, n_walkers, _ = samples.shape
    return n_steps * n_walkers / np.maximum(tau, 1.0)
