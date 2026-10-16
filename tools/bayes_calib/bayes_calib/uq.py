"""Propagate parameter uncertainty to the predictions the solver makes.

The solver returns ONE number per joint: Nf computed from fixed constants.
After calibration we hold thousands of plausible constant sets. Feeding the
same ΔW through each one gives a *distribution* of Nf per joint, so the
report can say "Nf = 2400 cycles (90 % credible 1900-3100)" and
"P(B10 life < required 1000 cycles) = 4 %", instead of a bare 2400.

Two sources of variation stay separate (the usual aleatory/epistemic split):

  * epistemic  - we don't know the true constants (width of the posterior).
                 Shrinks as you add test data.
  * aleatory   - identical joints still fail at different times (Weibull β).
                 Does NOT shrink with more data; it is physical scatter.

`eta` bands below show epistemic uncertainty only. `b10` folds in the
Weibull scatter for each draw, which is what a reliability target is
judged against.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .model import ln_characteristic_life


def b_life(eta: np.ndarray, beta: np.ndarray, p: float) -> np.ndarray:
    """Cycles by which fraction p of a population has failed:  η·(−ln(1−p))^(1/β)."""
    return eta * (-np.log1p(-p)) ** (1.0 / beta)


def predict_joints(
    theta_draws: np.ndarray,
    dw: np.ndarray,
    pad_d_mm: float,
    n_use: int = 500,
    seed: int = 0,
) -> dict[str, np.ndarray]:
    """Per-joint credible bands. Returns arrays of shape (n_joints,) per statistic."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(theta_draws), size=min(n_use, len(theta_draws)), replace=False)
    th = theta_draws[idx]
    d = np.full_like(dw, pad_d_mm, dtype=float)

    eta = np.exp(ln_characteristic_life(th, dw, d))          # (draws, joints)
    beta = np.exp(th[:, [4]])
    b10 = b_life(eta, beta, 0.10)
    b1 = b_life(eta, beta, 0.01)

    q = lambda a, p: np.percentile(a, p, axis=0)
    return {
        "eta_p05": q(eta, 5), "eta_p50": q(eta, 50), "eta_p95": q(eta, 95),
        "b10_p05": q(b10, 5), "b10_p50": q(b10, 50), "b10_p95": q(b10, 95),
        "b1_p05": q(b1, 5), "b1_p50": q(b1, 50), "b1_p95": q(b1, 95),
        "_b10_draws": b10,
    }


def propagate_solver_csv(
    theta_draws: np.ndarray,
    fatigue_csv: str | Path,
    pad_d_mm: float,
    out_csv: str | Path,
    target_cycles: float | None = None,
) -> int:
    """Read rpim_solver's `*_solder_fatigue.csv`, write per-joint UQ columns.

    Uses the `delta_W_mpa` column (the solver's creep-work output) and
    recomputes life under the posterior instead of the fixed constants.
    """
    with open(fatigue_csv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "delta_W_mpa" not in rows[0]:
        raise ValueError(f"{fatigue_csv}: no `delta_W_mpa` column (is this a solver fatigue CSV?)")

    dw = np.array([float(r["delta_W_mpa"]) for r in rows])
    pred = predict_joints(theta_draws, dw, pad_d_mm)

    extra = ["Nf_solver_cycles"] if "Nf_cycles" in rows[0] else []
    header = ["joint_id", "ref_designator", "pin_number", "delta_W_mpa", *extra,
              "eta_p05", "eta_p50", "eta_p95", "B10_p05", "B10_p50", "B10_p95",
              "B1_p05", "B1_p50", "B1_p95"]
    if target_cycles is not None:
        header.append(f"P(B10<{target_cycles:g})")

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i, r in enumerate(rows):
            line = [r.get("joint_id", i), r.get("ref_designator", ""), r.get("pin_number", ""),
                    f"{dw[i]:.6e}"]
            if extra:
                line.append(r["Nf_cycles"])
            line += [f"{pred[k][i]:.1f}" for k in
                     ("eta_p05", "eta_p50", "eta_p95", "b10_p05", "b10_p50", "b10_p95",
                      "b1_p05", "b1_p50", "b1_p95")]
            if target_cycles is not None:
                line.append(f"{np.mean(pred['_b10_draws'][:, i] < target_cycles):.4f}")
            w.writerow(line)
    return len(rows)
