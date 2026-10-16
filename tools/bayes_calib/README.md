# bayes_calib — Bayesian calibration & uncertainty quantification

A post-processing tool for ElectroniX. It never touches the solver. It reads
**fatigue test data** and **solver output (`ΔW`)**, and answers:

1. **What Darveaux constants (K1–K4) and Weibull shape (β) are consistent with
   my test data, and how sure are we?**  → *calibration*
2. **Given that uncertainty, what is the real range of life for each solder
   joint?**  → *uncertainty quantification (UQ)*

## What you get

| Output | Meaning |
|---|---|
| `summary.txt` / `posterior.json` | Median and 95 % credible interval for K1, K2, K3, K4, β, plus convergence diagnostics and warnings |
| `corner.png` | Posterior vs prior for every parameter, and the trade-offs between them |
| `prediction.png` | Life-vs-ΔW curve with a 90 % credible band, with your test data overlaid |
| `trace.png` | Sampler health check |
| `pcsim_override.txt` | Ready-to-paste `*MATERIAL_OVERRIDE` block (feeds the solver's existing override path) |
| `nf_uq.csv` (from `propagate`) | Per-joint η, B10 and B1 life with 5/50/95 % bands, and optionally `P(B10 < required cycles)` |

The headline change: instead of **"Nf = 2400 cycles"** you can report
**"Nf = 2400 (90 % credible 1900–3100); probability of missing the
1000-cycle requirement = 4 %"**.

## Quick start

```bash
cd tools/bayes_calib
pip install -r requirements.txt

python -m bayes_calib synth examples/synthetic_tests.csv      # make a demo campaign (known truth)
python -m bayes_calib calibrate examples/synthetic_tests.csv --out out
python -m bayes_calib propagate out/posterior_samples.npz \
        path/to/model_solder_fatigue.csv --pad-d 0.5 --target-cycles 1000 --out nf_uq.csv

python -m unittest discover -s tests -v                        # 7 tests, ~12 s
```

Your own data, one row per specimen:

```csv
group_id,delta_w_mpa,pad_d_mm,cycles,failed
TV1-A,0.015,0.50,4210,1
TV1-A,0.015,0.50,6000,0
```

`delta_w_mpa` is ΔW from an `rpim_solver` run of that test condition.
`failed=0` marks a specimen that survived until the test was stopped.

---

## Theory

### The problem

The solver predicts life from the Darveaux model:

```
N0    = K1 · ΔW^K2           cycles to crack initiation
da/dN = K3 · ΔW^K4           crack growth rate [mm/cycle]
Nf    = N0 + 0.65·d / (da/dN)   d = pad diameter
```

K1–K4 come from the literature for a given alloy and package, and the true
values for *your* board differ. A single set of constants hides that: the
solver prints one `Nf` with no indication of how far off it could be.

### Bayes' rule in one line

```
p(θ | data)  ∝  p(data | θ)  ×  p(θ)
 posterior      likelihood      prior
```

- **Prior `p(θ)`** — what we believed before the test: the library constants,
  with generous error bars (`model.py :: Prior`).
- **Likelihood `p(data | θ)`** — how probable the observed failure times are if
  the constants were `θ` (`model.py :: log_likelihood`).
- **Posterior `p(θ | data)`** — what we should believe now. It is a
  *distribution*, not a single number, and it narrows as evidence accumulates.

### The likelihood: Weibull with censoring

Identical joints don't fail at identical times; failure times follow a Weibull
distribution with characteristic life **η = Nf(θ)** and shape **β**.

- A **failed** specimen contributes the Weibull density `f(t)`.
- A **suspended** specimen (test stopped, still alive) contributes the
  survival probability `S(t) = exp(−(t/η)^β)`.

Dropping suspended specimens would bias life estimates *low*, because only the
early failures would be counted. Keeping them as censored data is correct.

### Sampling the posterior: MCMC

The posterior has no closed form, so we draw samples from it
(`sampler.py`) using the Goodman–Weare *affine-invariant ensemble sampler*:
many "walkers" explore the parameter space and each proposes moves along the
line to another walker. K1–K4 are strongly correlated (the posterior is a thin
ridge), and this sampler adapts to that shape without tuning. It is
implemented in ~60 lines of numpy so you can read it.

Once you hold ~15 000 samples, everything is an average:
the mean of K1 is the column mean, a 95 % interval is two percentiles, and the
uncertainty of any prediction is that prediction evaluated at every sample
(`uq.py`).

### Two kinds of uncertainty

| | Source | Shrinks with more data? | In the output |
|---|---|---|---|
| **Epistemic** | We don't know the true constants | **Yes** | width of the `η` band |
| **Aleatory** | Identical joints still fail at different times | **No** (physical) | Weibull β, folded into the `B10`/`B1` bands |

### Identifiability (read this before trusting individual constants)

The data constrain the *life curve*, not necessarily each constant. With the
solver's built-in SAC305 constants (K1=22.18, K3=0.2169), crack initiation is
91–96 % of life. Measured on synthetic data generated from exactly those
constants:

| param | true | posterior median | 95 % interval | shrinkage vs prior |
|---|---|---|---|---|
| K1 | 22.18 | 22.17 | 15.9 – 28.9 | 90 % |
| K2 | −1.26 | −1.267 | −1.349 – −1.191 | 92 % |
| **K3** | 0.2169 | 0.247 | **0.022 – 4.05** | **11 %** |
| K4 | 1.26 | 1.05 | 0.31 – 1.73 | 28 % |
| β | 5.4 | 4.80 | 4.18 – 5.46 | 86 % |

K3 and K4 describe the crack-growth phase, which is only ~5 % of life, so
total-life data are nearly blind to them. The tool prints a warning when this
happens. To pin them down you need data that separates the phases:
**crack-initiation observations** (in-situ resistance monitoring, dye-and-pry
at interim cycle counts) and a **spread of pad diameters**.

When the data are silent, the posterior simply equals the prior, which is the
honest answer.

### What the calibrated constants mean

`K1–K4` are fitted *together with this solver's ΔW*. If the solver's ΔW is
biased high by 20 %, calibration absorbs that into K1. The constants are
therefore valid **only with the same solver version, mesh density and
settings** that produced the ΔW in your data CSV. Do not carry them to
another tool.

---

## Validation of the tool itself

Because the truth is known for synthetic data, the tool can be tested honestly:

- **Coverage test** (24 independent synthetic campaigns): the 95 % interval
  contained the true value in 88–100 % of runs for every parameter, with no
  bias in β (mean posterior median 5.04 vs truth 5.0).
- **Sampler test**: recovers a correlated 2-D Gaussian's mean and covariance.
- **Likelihood test**: matches `scipy.stats.weibull_min` including censoring.
- **Forward model test**: reproduces `rpim_solver::darveaux_life` exactly.

## What this tool does NOT do (yet)

- It uses **no real data** out of the box. The bundled CSV is synthetic;
  results from it say nothing about SAC305 in the real world.
- ΔW is treated as exact. Solver error in ΔW is not propagated; adding a
  ΔW-bias parameter is the natural next step once real data exist.
- Initiation-only observations (`N0`) are not yet supported in the likelihood.
- Applies the single-pad-diameter assumption the solver itself makes.

## Files

```
bayes_calib/model.py       forward model, prior, likelihood, posterior
bayes_calib/sampler.py     ensemble MCMC + autocorrelation diagnostics
bayes_calib/calibrate.py   MAP search -> MCMC -> Posterior
bayes_calib/uq.py          propagate draws to per-joint life bands
bayes_calib/report.py      summaries, plots, .pcsim override
bayes_calib/data.py        CSV loader + synthetic campaign generator
tests/                     unit + end-to-end tests
```

## References

- Darveaux, R. (2002). *Effect of simulation methodology on solder joint crack
  growth correlation and fatigue life prediction*, J. Electron. Packag. 124(2).
- Goodman, J. & Weare, J. (2010). *Ensemble samplers with affine invariance*,
  Comm. App. Math. and Comp. Sci. 5(1).
- Gelman et al., *Bayesian Data Analysis*, 3rd ed. (censored data, §6; MCMC diagnostics, §11).
- JEDEC JESD22-A104 / IPC-9701A (thermal-cycling test and Weibull reporting).
