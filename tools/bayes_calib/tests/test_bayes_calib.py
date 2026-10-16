"""Run:  python -m unittest discover -s tests -v   (from tools/bayes_calib)"""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import weibull_min

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bayes_calib import calibrate, load_csv, make_synthetic, propagate_solver_csv
from bayes_calib.model import Dataset, Prior, ln_characteristic_life, log_likelihood
from bayes_calib.sampler import autocorr_time, ensemble_sample

TRUTH = np.array([np.log(22.18), -1.26, np.log(0.03), 1.10, np.log(5.0)])


class ForwardModel(unittest.TestCase):
    def test_matches_rust_darveaux_life(self):
        # Re-derive Nf with plain pow() exactly like rpim_solver::darveaux_life.
        k1, k2, k3, k4 = 22.18, -1.26, 0.2169, 1.26
        for dw, d in [(0.01, 0.5), (0.5, 0.3), (2.0, 0.76)]:
            n0 = k1 * dw ** k2
            rate = k3 * dw ** k4
            expected = n0 + 0.65 * d / rate
            theta = np.array([[np.log(k1), k2, np.log(k3), k4, 0.0]])
            got = np.exp(ln_characteristic_life(theta, np.array([dw]), np.array([d])))[0, 0]
            self.assertAlmostEqual(got / expected, 1.0, places=12)


class Likelihood(unittest.TestCase):
    def test_matches_scipy_weibull_with_censoring(self):
        rng = np.random.default_rng(3)
        dw = np.full(8, 0.03)
        d = np.full(8, 0.5)
        t = rng.uniform(500, 3000, 8)
        failed = np.array([1, 1, 1, 0, 1, 0, 1, 1], dtype=bool)
        data = Dataset(dw, d, t, failed, np.array(["g"] * 8))
        theta = TRUTH[None, :]

        eta = np.exp(ln_characteristic_life(theta, dw, d))[0]
        beta = np.exp(theta[0, 4])
        expected = (weibull_min.logpdf(t[failed], c=beta, scale=eta[failed]).sum()
                    + weibull_min.logsf(t[~failed], c=beta, scale=eta[~failed]).sum())
        self.assertAlmostEqual(log_likelihood(theta, data)[0], expected, places=9)


class Sampler(unittest.TestCase):
    def test_recovers_correlated_gaussian(self):
        mean = np.array([1.0, -2.0])
        cov = np.array([[1.0, 0.9], [0.9, 1.0]])   # the thin ridge MCMC finds hard
        prec = np.linalg.inv(cov)
        logp = lambda x: -0.5 * np.einsum("wi,ij,wj->w", x - mean, prec, x - mean)

        rng = np.random.default_rng(0)
        chain = ensemble_sample(logp, rng.standard_normal((20, 2)), 4000, seed=1)
        s = chain.samples[1000:].reshape(-1, 2)
        np.testing.assert_allclose(s.mean(0), mean, atol=0.06)
        np.testing.assert_allclose(np.cov(s.T), cov, atol=0.08)
        self.assertTrue(np.all(autocorr_time(chain.samples[1000:]) > 1.0))


class EndToEnd(unittest.TestCase):
    def test_posterior_predicts_true_life_within_credible_band(self):
        # Individual K's are partly degenerate, but the LIFE CURVE is what data
        # constrain. The true eta at every tested condition must fall inside the
        # posterior's 99 % band.
        data = make_synthetic(TRUTH, [0.015, 0.03, 0.06], [0.25, 0.8], 10, 6000.0, seed=7)
        post = calibrate(data, n_steps=2500, seed=2)
        conds = sorted({(a, b) for a, b in zip(data.dw, data.d_mm)})
        dw = np.array([c[0] for c in conds])
        d = np.array([c[1] for c in conds])
        eta_post = np.exp(ln_characteristic_life(post.theta[::5], dw, d))
        eta_true = np.exp(ln_characteristic_life(TRUTH[None, :], dw, d))[0]
        lo, hi = np.percentile(eta_post, [0.5, 99.5], axis=0)
        self.assertTrue(np.all((eta_true >= lo) & (eta_true <= hi)),
                        f"true eta outside band: {eta_true} vs [{lo}, {hi}]")

    def test_prior_is_respected_with_no_informative_data(self):
        # With almost no data the posterior must stay near the prior, not blow up.
        data = make_synthetic(TRUTH, [0.03], [0.5], 3, None, seed=1)
        data = Dataset(data.dw, data.d_mm, data.cycles, data.failed, data.group)
        post = calibrate(data, n_steps=2000, seed=3)
        prior = Prior()
        self.assertLess(abs(np.median(post.theta[:, 1]) - prior.mean[1]), 2 * prior.sd[1])
        self.assertTrue(np.all(post.theta[:, 1] < 0) and np.all(post.theta[:, 3] > 0))


class PropagateAndIO(unittest.TestCase):
    def test_propagate_writes_bands_and_probability(self):
        theta = TRUTH + 0.02 * np.random.default_rng(0).standard_normal((300, 5))
        with tempfile.TemporaryDirectory() as tmp:
            src, out = Path(tmp) / "f.csv", Path(tmp) / "uq.csv"
            with open(src, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["joint_id", "ref_designator", "pin_number", "delta_W_mpa", "Nf_cycles"])
                w.writerow([0, "U1", "A1", 0.02, 3000])
                w.writerow([1, "U1", "A2", 0.06, 900])
            n = propagate_solver_csv(theta, src, 0.5, out, target_cycles=500)
            self.assertEqual(n, 2)
            with open(out, newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertLess(float(rows[0]["eta_p05"]), float(rows[0]["eta_p95"]))
            self.assertGreater(float(rows[0]["eta_p50"]), float(rows[1]["eta_p50"]))  # higher dW, shorter life
            self.assertIn("P(B10<500)", rows[0])

    def test_load_csv_rejects_bad_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.csv"
            bad.write_text("delta_w_mpa,pad_d_mm\n0.1,0.5\n")
            with self.assertRaisesRegex(ValueError, "missing column"):
                load_csv(bad)


if __name__ == "__main__":
    unittest.main()
