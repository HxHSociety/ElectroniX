"""Command line:  python -m bayes_calib <synth|calibrate|propagate> ..."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from .calibrate import calibrate
from .data import load_csv, make_synthetic, save_csv
from .report import text_summary, write_all
from .uq import propagate_solver_csv

# Illustrative synthetic truth. NOT a material claim: chosen so that crack
# growth is a visible share of life, which makes K3 learnable in the demo.
SYNTH_TRUTH = np.array([np.log(22.18), -1.26, np.log(0.03), 1.10, np.log(5.0)])


def cmd_synth(a: argparse.Namespace) -> int:
    data = make_synthetic(
        SYNTH_TRUTH,
        dw_levels=[0.015, 0.025, 0.04, 0.07],
        d_levels=[0.25, 0.5, 0.8],
        n_per_group=a.n,
        test_end_cycles=a.test_end,
        seed=a.seed,
    )
    save_csv(data, a.out)
    print(f"wrote {a.out}: {len(data)} specimens, {int((~data.failed).sum())} censored")
    print("truth: K1=22.18 K2=-1.26 K3=0.03 K4=1.10 beta=5.0")
    return 0


def cmd_calibrate(a: argparse.Namespace) -> int:
    data = load_csv(a.data)
    t0 = time.time()
    print(f"calibrating on {len(data)} specimens ({a.walkers} walkers x {a.steps} steps)...")
    post = calibrate(data, n_walkers=a.walkers, n_steps=a.steps, seed=a.seed)
    print(f"sampling done in {time.time() - t0:.1f}s\n")
    print(text_summary(post, data))
    files = write_all(post, data, a.out, material_id=a.material_id)
    print("\nwrote:")
    for f in files:
        print(f"  {f}")
    return 0


def cmd_propagate(a: argparse.Namespace) -> int:
    theta = np.load(a.samples)["theta"]
    n = propagate_solver_csv(theta, a.fatigue_csv, a.pad_d, a.out, a.target_cycles)
    print(f"wrote {a.out}: {n} joints with 90 % credible bands")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="bayes_calib", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("synth", help="generate a synthetic test campaign (known truth)")
    s.add_argument("--out", type=Path, default=Path("examples/synthetic_tests.csv"))
    s.add_argument("--n", type=int, default=10, help="specimens per condition")
    s.add_argument("--test-end", type=float, default=6000.0, help="suspend tests at this many cycles")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_synth)

    c = sub.add_parser("calibrate", help="infer Darveaux constants from test data")
    c.add_argument("data", type=Path, help="CSV: group_id,delta_w_mpa,pad_d_mm,cycles,failed")
    c.add_argument("--out", type=Path, default=Path("out"))
    c.add_argument("--steps", type=int, default=6000)
    c.add_argument("--walkers", type=int, default=40)
    c.add_argument("--seed", type=int, default=1)
    c.add_argument("--material-id", type=int, default=4, help="ID for *MATERIAL_OVERRIDE (4=SAC305)")
    c.set_defaults(fn=cmd_calibrate)

    u = sub.add_parser("propagate", help="credible bands on a solver fatigue CSV")
    u.add_argument("samples", type=Path, help="posterior_samples.npz from `calibrate`")
    u.add_argument("fatigue_csv", type=Path, help="rpim_solver *_solder_fatigue.csv")
    u.add_argument("--pad-d", type=float, required=True, help="pad diameter used in that solve [mm]")
    u.add_argument("--target-cycles", type=float, default=None,
                   help="also report P(B10 life < target)")
    u.add_argument("--out", type=Path, default=Path("nf_uq.csv"))
    u.set_defaults(fn=cmd_propagate)

    args = p.parse_args(argv)
    try:
        return args.fn(args)
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
