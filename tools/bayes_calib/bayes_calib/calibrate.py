"""End-to-end calibration: find a good start, run MCMC, return posterior samples."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .model import Dataset, Prior, make_log_posterior, to_natural
from .sampler import Chain, autocorr_time, effective_sample_size, ensemble_sample


@dataclass
class Posterior:
    theta: np.ndarray          # (n_draws, 5)  [ln K1, K2, ln K3, K4, ln β]
    log_prob: np.ndarray       # (n_draws,)
    theta_map: np.ndarray      # (5,)  highest-posterior draw (a coherent point)
    chain: Chain               # full chain, for trace plots
    burn: int
    tau: np.ndarray            # autocorrelation time per parameter
    ess: np.ndarray            # effective sample size per parameter
    prior: Prior

    @property
    def natural(self) -> np.ndarray:
        """Draws as [K1, K2, K3, K4, β]."""
        return to_natural(self.theta)

    def acceptance_bad(self) -> bool:
        """Stretch-move samplers are healthy at roughly 20-50 % acceptance."""
        return not (0.1 <= self.chain.acceptance <= 0.7)


def find_map(log_post, prior: Prior, rng: np.random.Generator, n_starts: int = 12) -> np.ndarray:
    """Maximum-a-posteriori start point via multi-start Nelder-Mead.

    MCMC converges far faster when the walkers start near the posterior mass
    rather than scattered across the prior.
    """
    best, best_val = None, np.inf
    for k in range(n_starts):
        x0 = prior.mean if k == 0 else prior.mean + prior.sd * 0.5 * rng.standard_normal(5)
        if not np.isfinite(log_post(x0[None, :])[0]):
            continue
        res = minimize(lambda x: -log_post(x[None, :])[0], x0,
                       method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000})
        if np.isfinite(res.fun) and res.fun < best_val:
            best, best_val = res.x, res.fun
    if best is None:
        raise RuntimeError("could not find any finite-posterior starting point")
    return best


def calibrate(
    data: Dataset,
    prior: Prior | None = None,
    n_walkers: int = 40,
    n_steps: int = 6000,
    burn_frac: float = 0.4,
    thin: int = 10,
    seed: int = 1,
) -> Posterior:
    prior = prior or Prior()
    rng = np.random.default_rng(seed)
    log_post = make_log_posterior(data, prior)

    x_map = find_map(log_post, prior, rng)

    # Tight ball of walkers around the MAP; reject any non-finite starts.
    p0 = x_map + 1e-3 * rng.standard_normal((n_walkers, 5))
    for i in range(n_walkers):
        while not np.isfinite(log_post(p0[i:i + 1])[0]):
            p0[i] = x_map + 1e-4 * rng.standard_normal(5)

    chain = ensemble_sample(log_post, p0, n_steps, seed=seed)

    burn = int(burn_frac * n_steps)
    kept = chain.samples[burn:]
    tau = autocorr_time(kept)
    ess = effective_sample_size(kept, tau)

    flat = kept[::thin].reshape(-1, 5)
    flat_lp = chain.log_prob[burn:][::thin].reshape(-1)
    return Posterior(
        theta=flat, log_prob=flat_lp, theta_map=flat[int(np.argmax(flat_lp))],
        chain=chain, burn=burn, tau=tau, ess=ess, prior=prior,
    )
