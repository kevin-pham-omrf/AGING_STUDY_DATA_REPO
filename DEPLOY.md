# Deploying the app

The app runs as a single container: `python:3.13-slim`, a non-root user, hash-checked
dependencies, and the app code and data files baked in. Expect an image of roughly 0.8 GB
(about 630 MB of Python packages and 106 MB of app and data, plus the base image).

Before this is reachable from the internet, its hosting has to be approved through Research
Computing / IT Security; don't run it from the HPC systems.

## Run it

```bash
docker compose up -d --build
docker compose ps        # STATUS shows "healthy" once the app has loaded its data (~15 s)
```

The app is then at <http://127.0.0.1:8000>. `compose.yaml` publishes it on localhost only and
runs it with a read-only filesystem, no Linux capabilities, `no-new-privileges`, and a
1.5 GB memory cap. Without Compose, the equivalent is:

```bash
docker build -t aging-study-app .
docker run -d --name aging-study-app --restart unless-stopped -p 127.0.0.1:8000:8000 \
  --read-only --tmpfs /tmp --cap-drop ALL --security-opt no-new-privileges \
  --memory 1536m --cpus 1 --pids-limit 128 aging-study-app
```

## Putting it on the internet

Put a reverse proxy in front that terminates TLS. It must pass websocket upgrades (the app
talks to the browser over a websocket) and allow long-lived idle connections. An nginx example:

```nginx
limit_req_zone $binary_remote_addr zone=aging:10m rate=10r/s;

map $http_upgrade $connection_upgrade {
    default upgrade;
    ''      close;
}

server {
    listen 443 ssl;
    server_name aging.example.org;                  # your hostname
    ssl_certificate     /etc/ssl/certs/aging.pem;   # your certificate
    ssl_certificate_key /etc/ssl/private/aging.key;

    client_max_body_size 1m;   # the app takes no uploads

    location / {
        limit_req zone=aging burst=40 nodelay;      # a page load fetches ~20 files
        proxy_pass http://127.0.0.1:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection $connection_upgrade;
        proxy_set_header Host $host;
        proxy_read_timeout 1h;                      # visitors may sit on a plot for a while
        proxy_buffering off;
    }
}
```

`limit_req` limits HTTP requests (page loads, downloads, new websocket connections), not the
messages inside an open websocket. Those are limited by how fast the app answers them (about
0.3 s of server time per gene change), by `--ws-max-size` (1 MiB per message), and by
`--limit-concurrency` (200 open connections per container, after which new requests get a 503).
Both flags are set in the Dockerfile's `CMD`.

## Sizing and scaling

- About 0.6 GB of memory per container, ~15 s to start (all data is loaded into memory), and
  one CPU. Each container runs one app process.
- A gene change costs ~0.3 s of server time with every plot on screen; in a browser only the
  visible tab's plots are drawn, so usually less. One container comfortably serves a lab-sized
  audience. While one visitor's request is being handled, others wait, so a burst of heavy use
  shows up as slower responses for everyone on that container.
- To scale, run more containers behind the proxy with **sticky sessions** (nginx `ip_hash`, or
  a cookie-based option). A visitor's websocket and their CSV downloads must reach the same
  process, so don't use `uvicorn --workers` or round-robin load balancing.

## Updating

- **Base image security updates:** rebuild regularly with `docker compose build --pull`, then
  `docker compose up -d`. To make builds reproducible, pin the base image by digest in the
  Dockerfile (`FROM python:3.13-slim@sha256:...`) and update the digest deliberately.
- **Python dependencies:** edit `requirements.in`, then regenerate the hash-locked
  `requirements.txt` with Python 3.13 and audit it:

  ```bash
  uv pip compile requirements.in --python-version 3.13 --generate-hashes -o requirements.txt
  pip-audit -r requirements.txt
  ```

- **Data:** replace the parquet files and rebuild. The app reads
  `ALL_TPM_LABELED_FILTERED.parquet`, `ALL_GENE_BODY_PER_SAMPLE_T_v2.parquet`,
  `ALL_TSS_PER_SAMPLE_T_v2.parquet` and `DETECTED_GENES.csv`.
- **Testing a build:** the test suite can check a running container against any git revision
  (see `tests/README.md`):

  ```bash
  CAND_URL=http://127.0.0.1:8000 CAND_PYTHON=/path/to/venv/bin/python tests/run.sh
  ```

## Known follow-ups

- Shiny 1.8 logs a deprecation warning for `render.download` (three per visitor session). It
  works, but should move to `render.download_button` before upgrading Shiny.
- The footer's link to the study is empty.
- Promoter methylation values are clipped at 10% in the box plots, gene body values at 0%.
