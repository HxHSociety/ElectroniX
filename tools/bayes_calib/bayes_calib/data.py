"""Reading real test data and generating synthetic test data.

Real data CSV (one row per specimen):

    group_id,delta_w_mpa,pad_d_mm,cycles,failed
    TV1-A,0.015,0.50,4210,1
    TV1-A,0.015,0.50,5100,1
    TV1-A,0.015,0.50,6000,0        <- test stopped at 6000 cycles, never failed

* delta_w_mpa : ΔW from an rpim_solver run of THAT test condition
                (e.g. the average ΔW of the failed joint's neighbourhood, or
                 the `delta_W_mpa` column of the solver's fatigue CSV).
* pad_d_mm    : pad diameter used in that solver run.
* cycles      : cycles-to-failure, or cycles at which the test was suspended.
* failed      : 1 = failure observed, 0 = right-censored (survived).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .model import Dataset, ln_characteristic_life


def load_csv(path: str | Path) -> Dataset:
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        required = {"delta_w_mpa", "pad_d_mm", "cycles", "failed"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for i, r in enumerate(reader, start=2):
            try:
                rows.append((
                    r.get("group_id") or f"row{i}",
                    float(r["delta_w_mpa"]), float(r["pad_d_mm"]),
                    float(r["cycles"]), int(r["failed"]),
                ))
            except ValueError as e:
                raise ValueError(f"{path}:{i}: {e}") from e

    if len(rows) < 5:
        raise ValueError("need at least 5 specimens to calibrate anything")
    group, dw, d, t, f = zip(*rows)
    dw, d, t = (np.array(v, dtype=float) for v in (dw, d, t))
    if np.any(dw <= 0) or np.any(d <= 0) or np.any(t <= 0):
        raise ValueError("delta_w_mpa, pad_d_mm and cycles must all be > 0")
    failed = np.array(f, dtype=bool)
    if not failed.any():
        raise ValueError("every specimen is censored; no failures to learn from")
    return Dataset(dw=dw, d_mm=d, cycles=t, failed=failed, group=np.array(group))


def make_synthetic(
    theta_true: np.ndarray,
    dw_levels: list[float],
    d_levels: list[float],
    n_per_group: int = 10,
    test_end_cycles: float | None = 6000.0,
    seed: int = 0,
) -> Dataset:
    """Simulate a fatigue test campaign from KNOWN parameters.

    Weibull sampling by inverse CDF:  t = eta · (−ln U)^(1/beta),  U ~ Uniform(0,1).
    Specimens still alive at `test_end_cycles` are recorded as censored, exactly
    like a real test that is stopped for schedule reasons.

    Because the truth is known, we can check that calibration recovers it.
    That is the only honest way to test a calibration tool before real data exist.
    """
    rng = np.random.default_rng(seed)
    beta = float(np.exp(theta_true[4]))
    dw, d, t, failed, group = [], [], [], [], []
    for dwi in dw_levels:
        for di in d_levels:
            ln_eta = ln_characteristic_life(
                np.asarray(theta_true)[None, :], np.array([dwi]), np.array([di]))[0, 0]
            eta = np.exp(ln_eta)
            life = eta * (-np.log(rng.random(n_per_group))) ** (1.0 / beta)
            for ti in life:
                dead = test_end_cycles is None or ti <= test_end_cycles
                dw.append(dwi); d.append(di)
                t.append(ti if dead else test_end_cycles)
                failed.append(dead)
                group.append(f"dW{dwi:g}_d{di:g}")
    return Dataset(np.array(dw), np.array(d), np.array(t),
                   np.array(failed, dtype=bool), np.array(group))


def save_csv(data: Dataset, path: str | Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["group_id", "delta_w_mpa", "pad_d_mm", "cycles", "failed"])
        for g, dw, d, t, f in zip(data.group, data.dw, data.d_mm, data.cycles, data.failed):
            w.writerow([g, f"{dw:.6g}", f"{d:.6g}", f"{t:.1f}", int(f)])
