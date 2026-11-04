# Aging study data browser (Shiny for Python). See DEPLOY.md.
FROM python:3.14-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Unprivileged user to run the app; it only needs to read its files
RUN useradd --system --uid 10001 --no-create-home --shell /usr/sbin/nologin app

WORKDIR /app

# Every package is pinned and hash-checked (requirements.txt is generated from requirements.in)
COPY requirements.txt .
RUN pip install --require-hashes -r requirements.txt

# App code and data, owned by root and read-only to the app user
COPY app.py DETECTED_GENES.csv ./
COPY www/ www/
COPY ALL_TPM_LABELED_FILTERED.parquet ALL_GENE_BODY_PER_SAMPLE_T_v2.parquet ALL_TSS_PER_SAMPLE_T_v2.parquet ./

USER app
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/', timeout=4)"]

# A single worker: a Shiny session (its websocket and its download requests) must stay on
# the process that created it. Scale by running more containers behind a proxy with sticky
# sessions (see DEPLOY.md).
#   --ws-max-size: browser-to-server messages are small input updates; cap them at 1 MiB
#   --limit-concurrency: answer 503 instead of queueing once this many connections are open
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1", \
     "--ws-max-size", "1048576", "--limit-concurrency", "200", "--no-server-header"]
