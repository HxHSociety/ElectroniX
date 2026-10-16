FROM python:3.12-slim

# Headless plotting; never buffer logs; no .pyc clutter.
ENV MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Dependencies first so this layer is cached until requirements change.
COPY requirements.txt requirements-service.txt ./
RUN pip install -r requirements.txt -r requirements-service.txt

COPY bayes_calib ./bayes_calib
COPY service ./service

# Run as an unprivileged user.
RUN useradd --create-home --uid 10001 app
USER app

# Cloud Run (and most platforms) inject $PORT; default to 8080 locally.
# One worker: each request is CPU-heavy, scale by adding instances instead.
CMD ["sh", "-c", "exec uvicorn service.app:app --host 0.0.0.0 --port ${PORT:-8080} --workers 1"]
