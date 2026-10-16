"""HTTP API around bayes_calib.

    POST /v1/calibrate?steps=6000&walkers=40&material_id=4&include_plots=false
         body: the test-data CSV, sent raw (Content-Type: text/csv)
    GET  /healthz

Design notes
------------
* The request is synchronous. A typical run takes 5-20 s, which fits inside a
  serverless request timeout. A job queue is the upgrade path for longer runs.
* MCMC is CPU work that any caller could abuse, so the API enforces limits
  BEFORE running: body size, row count, and a total-cost cap
  (steps x walkers x rows, calibrated on this tool's ~2.6 million units/s).
* Auth is a single shared key from the API_KEY env var. If API_KEY is unset
  the service runs open, which is only appropriate for local development.
"""

from __future__ import annotations

import base64
import os
import secrets
import tempfile
from pathlib import Path

from fastapi import Depends, FastAPI, Header, HTTPException, Query, Request
from fastapi.concurrency import run_in_threadpool

from bayes_calib import calibrate, load_csv
from bayes_calib.report import pcsim_override, summarize, text_summary, write_all

MAX_BODY_BYTES = 1_000_000        # 1 MB of CSV is ~20 000 specimens
MAX_ROWS = 5_000
MAX_COST = 150_000_000            # steps x walkers x rows  (~60 s of compute)
MAX_STEPS = 20_000
MAX_WALKERS = 80

app = FastAPI(title="ElectroniX Bayesian Calibration", version="0.1.0")


def require_key(x_api_key: str | None = Header(default=None)) -> None:
    expected = os.environ.get("API_KEY")
    if not expected:
        return  # open: local development only
    if not x_api_key or not secrets.compare_digest(x_api_key, expected):
        raise HTTPException(status_code=401, detail="missing or invalid X-API-Key")


@app.get("/healthz")
def healthz() -> dict:
    return {"status": "ok"}


@app.post("/v1/calibrate", dependencies=[Depends(require_key)])
async def calibrate_endpoint(
    request: Request,
    steps: int = Query(6000, ge=500, le=MAX_STEPS),
    walkers: int = Query(40, ge=10, le=MAX_WALKERS),
    material_id: int = Query(4, ge=1, le=999),
    seed: int = Query(1),
    include_plots: bool = Query(False),
) -> dict:
    if walkers % 2:
        raise HTTPException(422, "walkers must be even")

    body = await request.body()
    if len(body) > MAX_BODY_BYTES:
        raise HTTPException(413, f"CSV larger than {MAX_BODY_BYTES} bytes")
    if not body.strip():
        raise HTTPException(422, "empty body: send the test-data CSV")

    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "tests.csv"
        src.write_bytes(body)
        try:
            data = load_csv(src)
        except ValueError as e:
            raise HTTPException(422, f"bad CSV: {e}") from e
        if len(data) > MAX_ROWS:
            raise HTTPException(413, f"more than {MAX_ROWS} specimens")

        cost = steps * walkers * len(data)
        if cost > MAX_COST:
            raise HTTPException(
                413, f"job too large (cost {cost:.2e} > {MAX_COST:.2e}); reduce steps or walkers")

        post = await run_in_threadpool(
            calibrate, data, None, walkers, steps, 0.4, 10, seed)

        result: dict = {
            "n_specimens": len(data),
            "n_failed": int(data.failed.sum()),
            "acceptance": post.chain.acceptance,
            "summary": summarize(post),
            "text": text_summary(post, data),
            "pcsim_override": pcsim_override(post, material_id),
        }
        if include_plots:
            await run_in_threadpool(write_all, post, data, Path(tmp) / "out", material_id)
            result["plots_png_base64"] = {
                name: base64.b64encode((Path(tmp) / "out" / f"{name}.png").read_bytes()).decode()
                for name in ("corner", "trace", "prediction")
            }
    return result
