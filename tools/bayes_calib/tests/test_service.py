"""API tests. Needs: pip install fastapi uvicorn httpx"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fastapi.testclient import TestClient  # noqa: E402

from service.app import app  # noqa: E402

CSV = (ROOT / "examples" / "synthetic_tests.csv").read_bytes()
client = TestClient(app)


class Service(unittest.TestCase):
    def test_health(self):
        self.assertEqual(client.get("/healthz").json(), {"status": "ok"})

    def test_calibrate_returns_summary_and_override(self):
        r = client.post("/v1/calibrate?steps=1500&walkers=30", content=CSV,
                        headers={"Content-Type": "text/csv"})
        self.assertEqual(r.status_code, 200, r.text)
        j = r.json()
        self.assertEqual(j["n_specimens"], 120)
        self.assertIn("K1", j["summary"])
        self.assertIn("*MATERIAL_OVERRIDE, ID=4", j["pcsim_override"])
        self.assertNotIn("plots_png_base64", j)

    def test_plots_are_optional(self):
        r = client.post("/v1/calibrate?steps=1000&walkers=20&include_plots=true", content=CSV)
        self.assertEqual(r.status_code, 200, r.text)
        png = r.json()["plots_png_base64"]["prediction"]
        self.assertTrue(png.startswith("iVBOR"))   # base64 of the PNG magic bytes

    def test_rejects_bad_csv(self):
        r = client.post("/v1/calibrate", content=b"delta_w_mpa,pad_d_mm\n0.1,0.5\n")
        self.assertEqual(r.status_code, 422)
        self.assertIn("missing column", r.text)

    def test_rejects_empty_body(self):
        self.assertEqual(client.post("/v1/calibrate", content=b"").status_code, 422)

    def test_rejects_oversized_job_before_running(self):
        r = client.post("/v1/calibrate?steps=20000&walkers=80", content=CSV)
        self.assertEqual(r.status_code, 413)
        self.assertIn("too large", r.text)

    def test_rejects_out_of_range_parameters(self):
        self.assertEqual(client.post("/v1/calibrate?steps=10", content=CSV).status_code, 422)
        self.assertEqual(client.post("/v1/calibrate?walkers=41", content=CSV).status_code, 422)

    def test_api_key_enforced_when_configured(self):
        with mock.patch.dict(os.environ, {"API_KEY": "s3cret"}):
            url = "/v1/calibrate?steps=1000&walkers=20"
            self.assertEqual(client.post(url, content=CSV).status_code, 401)
            self.assertEqual(client.post(url, content=CSV, headers={"X-API-Key": "nope"}).status_code, 401)
            self.assertEqual(client.post(url, content=CSV, headers={"X-API-Key": "s3cret"}).status_code, 200)
            self.assertEqual(client.get("/healthz").status_code, 200)   # health stays open


if __name__ == "__main__":
    unittest.main()
