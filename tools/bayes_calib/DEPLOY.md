# Deploying the calibration service

The API lives in `service/app.py` and is packaged by `Dockerfile`.
All commands run from `tools/bayes_calib/`.

> Status: the API is covered by tests and was exercised against a live local
> server. The **Docker build and the cloud deploy below have not been run yet**
> (Docker and gcloud were not installed on the development machine). Expect to
> fix small issues on the first attempt.

## 1. Run it locally (no Docker)

```bash
pip install -r requirements.txt -r requirements-service.txt
API_KEY=dev python -m uvicorn service.app:app --port 8080

curl -X POST -H "X-API-Key: dev" -H "Content-Type: text/csv" \
     --data-binary @examples/synthetic_tests.csv \
     "http://localhost:8080/v1/calibrate?steps=6000"
```

Add `&include_plots=true` to also get `corner`, `trace` and `prediction` PNGs
as base64 in the response. Interactive docs: `http://localhost:8080/docs`.

## 2. Run it in Docker

```bash
docker build -t bayes-calib .
docker run --rm -p 8080:8080 -e API_KEY=dev bayes-calib
```

## 3. Deploy to Google Cloud Run

Cloud Run runs the container on demand and scales to zero when idle, so you
pay only while requests are running.

```bash
# one-time setup
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
gcloud services enable run.googleapis.com cloudbuild.googleapis.com \
                       artifactregistry.googleapis.com secretmanager.googleapis.com

# store the API key as a secret (never put it in the command line or the image)
openssl rand -hex 24 | tr -d '\n' | gcloud secrets create bayes-calib-key --data-file=-
# grant the Cloud Run runtime service account access to it:
gcloud secrets add-iam-policy-binding bayes-calib-key \
  --member="serviceAccount:$(gcloud projects describe YOUR_PROJECT_ID --format='value(projectNumber)')-compute@developer.gserviceaccount.com" \
  --role=roles/secretmanager.secretAccessor

# build from the Dockerfile and deploy
gcloud run deploy bayes-calib \
  --source . \
  --region us-central1 \
  --allow-unauthenticated \
  --set-secrets API_KEY=bayes-calib-key:latest \
  --cpu 2 --memory 1Gi \
  --concurrency 1 \
  --timeout 120 \
  --max-instances 3
```

Why these flags:

| Flag | Reason |
|---|---|
| `--allow-unauthenticated` | Cloud Run's own IAM is bypassed; the app enforces `X-API-Key` instead, so the desktop app can call it without Google credentials. |
| `--concurrency 1` | One calibration per instance: each request saturates a CPU core. Extra requests start extra instances. |
| `--max-instances 3` | Hard cost ceiling. Raise it deliberately. |
| `--timeout 120` | The API's own cost cap keeps jobs near 60 s. |

Test the deployment:

```bash
URL=$(gcloud run services describe bayes-calib --region us-central1 --format='value(status.url)')
KEY=$(gcloud secrets versions access latest --secret=bayes-calib-key)
curl "$URL/healthz"
curl -X POST -H "X-API-Key: $KEY" -H "Content-Type: text/csv" \
     --data-binary @examples/synthetic_tests.csv "$URL/v1/calibrate"
```

## AWS equivalent

Build and push the same image to ECR, then run it on **App Runner** (closest
to Cloud Run) or ECS Fargate. Put `API_KEY` in Secrets Manager and map it to
an environment variable. Set App Runner's max concurrency to 1.

## Before you put real customer data through it

- Test data and boards are customer IP. The service stores **nothing**: the CSV
  is processed in memory and a temp directory that is deleted on return. Keep it
  that way unless you add storage deliberately, and then add encryption,
  retention limits and per-tenant isolation.
- A single shared API key means no per-user quotas or revocation. For a real
  product, move to per-user keys or an identity provider (Clerk, Auth0).
- Browsers/Tauri calling this need the service URL added to `connect-src` in
  `src-tauri/tauri.conf.json`, which currently blocks external requests on
  purpose.
- Put a billing alert on the project before exposing the URL.
