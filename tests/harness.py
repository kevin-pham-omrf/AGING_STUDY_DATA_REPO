"""Drive a reference and a candidate version of the app over Shiny's websocket protocol.

Run via tests/run.sh, which starts both apps and sets REF_URL / CAND_URL.

    python tests/harness.py downloads|inputs|darkmode|perf

Exits non-zero if a check fails (perf only reports).
"""
import asyncio
import json
import os
import random
import sys
import time
import urllib.request

import pandas as pd
import pyarrow.parquet as pq
import websockets

REF_URL = os.environ.get("REF_URL", "http://127.0.0.1:8101")
CAND_URL = os.environ.get("CAND_URL", "http://127.0.0.1:8102")
APPS = [("reference", REF_URL), ("candidate", CAND_URL)]

WIDGETS = ["expression_plot", "tss_plot", "gene_body_plot", "tss_corr_plot", "gene_corr_plot"]
ALL_SELECTED = {"age": ["Young", "Adult", "Old"], "sex": ["Female", "Male"],
                "line": ["Microglia", "Neurons", "Astrocytes"]}
BASE = {"gene": "Cx3cr1", "filter": "7", **ALL_SELECTED, "toggle_dark:shiny.action": 0,
        **{f".clientdata_output_{w}_hidden": False for w in WIDGETS}}
SANITIZED = "An error has occurred. Check your logs or contact the app author for clarification."


class Client:
    """A minimal Shiny browser session: sets inputs, collects outputs, fetches downloads."""
    quiet = 1.5     # seconds of silence after an "idle" message that ends a round trip
    timeout = 60    # seconds to wait for the first "idle"

    async def open(self, url, **inputs):
        self.url = url
        ws_url = url.replace("http://", "ws://").replace("https://", "wss://")
        self.ws = await websockets.connect(f"{ws_url}/websocket/", max_size=None)
        self.sid = json.loads(await self.ws.recv())["config"]["sessionId"]
        return await self.send("init", {**BASE, **inputs})

    async def send(self, method, data):
        """Send inputs; collect output values, errors and custom messages until idle.

        "secs" is the time from sending to the last "idle" message (server work), not
        including the quiet-period wait.
        """
        t = time.perf_counter()
        await self.ws.send(json.dumps({"method": method, "data": data}))
        errors, values, custom, last_idle = {}, {}, [], None
        # Shiny may flush several times per update; read until quiet after an idle
        while True:
            try:
                msg = json.loads(await asyncio.wait_for(self.ws.recv(), self.quiet if last_idle else self.timeout))
            except asyncio.TimeoutError:
                break
            errors.update(msg.get("errors") or {})
            values.update(msg.get("values") or {})
            if "custom" in msg:
                custom.append(msg["custom"])
            if msg.get("busy") == "idle":
                last_idle = time.perf_counter()
        return {"errors": errors, "values": values, "custom": custom,
                "secs": (last_idle or time.perf_counter()) - t}

    def download(self, name):
        url = f"{self.url}/session/{self.sid}/download/{name}?w="
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            return f"<FAILED: {type(e).__name__}>"

    async def close(self):
        await self.ws.close()


def genes_with_all_data(n, seed):
    """`n` random detected genes that have expression and both kinds of modification data."""
    common = set(pd.read_csv("DETECTED_GENES.csv", header=None)[0])
    for f in ["ALL_TPM_LABELED_FILTERED.parquet", "ALL_GENE_BODY_PER_SAMPLE_T_v2.parquet",
              "ALL_TSS_PER_SAMPLE_T_v2.parquet"]:
        common &= set(pq.read_schema(f).names)
    random.seed(seed)
    return random.sample(sorted(common), n)


def short(errors):
    return {k: v["message"][:90] + (f"... ({len(v['message'])} chars)" if len(v["message"]) > 90 else "")
            for k, v in errors.items()}


async def check_downloads():
    """The three CSV downloads must match exactly across genes, groupings and filters."""
    genes = ["Cx3cr1"] + genes_with_all_data(2, seed=1)
    selections = [{"age": ["Old"], "sex": ["Male"], "line": ["Neurons", "Microglia"]},
                  {"age": ["Young", "Old"], "sex": ["Female"], "line": ["Astrocytes"]}]
    ref, cand = Client(), Client()
    ref.quiet = cand.quiet = 0.3
    ref.timeout = cand.timeout = 2  # plots hidden: downloads don't depend on them
    hidden = {f".clientdata_output_{w}_hidden": True for w in WIDGETS}
    await ref.open(REF_URL, **hidden)
    await cand.open(CAND_URL, **hidden)
    checked = mismatched = 0
    for g in genes:
        for flt in "1234567":
            for sel in selections:
                update = {"gene": g, "filter": flt, **ALL_SELECTED, **sel}
                await ref.send("update", update)
                await cand.send("update", update)
                for d in ["download_expr", "download_gene", "download_tss"]:
                    a, b = ref.download(d), cand.download(d)
                    checked += 1
                    if a != b or a.startswith("<FAILED"):
                        mismatched += 1
                        print(f"  MISMATCH {g} grouping {flt} {sel} {d}: {len(a)} vs {len(b)} chars")
    print(f"downloads compared: {checked}, mismatched: {mismatched}  (genes: {genes})")
    return mismatched == 0


async def check_inputs():
    """Forged inputs and genes without data must never show raw exception text."""
    genes = pd.read_csv("DETECTED_GENES.csv", header=None)[0].tolist()
    tss_cols = set(pq.read_schema("ALL_TSS_PER_SAMPLE_T_v2.parquet").names)
    expr_cols = set(pq.read_schema("ALL_TPM_LABELED_FILTERED.parquet").names)
    missing = next(g for g in genes if g not in tss_cols and g in expr_cols)
    good = {"gene": "Cx3cr1", "filter": "7", **ALL_SELECTED}
    cases = [("forged gene 'not_a_gene'", {**good, "gene": "not_a_gene"}),
             ("forged gene as list", {**good, "gene": ["a", "b"]}),
             ("forged filter '8'", {**good, "filter": "8"}),
             ("forged age as string", {**good, "age": "Young"}),
             (f"real gene without promoter data ({missing})", {**good, "gene": missing})]
    ok = True
    for name, url in APPS:
        print(name)
        c = Client()
        await c.open(url)
        for label, update in cases:
            r = await c.send("update", update)
            print(f"  {label}: errors={short(r['errors'])}")
            leaked = [m["message"] for m in r["errors"].values()
                      if m["message"] != SANITIZED and not m["message"].startswith("No ")]
            if name == "candidate" and leaked:
                ok = False
                print(f"    FAIL: raw error text shown: {leaked[0][:120]!r}")
        await c.close()
    return ok


async def check_darkmode():
    """One visitor's dark-mode toggle must not affect another visitor's."""
    ok = True
    for name, url in APPS:
        a, b = Client(), Client()
        await a.open(url)
        await b.open(url)
        modes = []
        for c in (a, b):
            r = await c.send("update", {"toggle_dark:shiny.action": 1})
            modes.append([m["bslib.toggle-dark-mode"]["value"] for m in r["custom"] if "bslib.toggle-dark-mode" in m])
        print(f"{name}: visitor A first click -> {modes[0]}, visitor B first click -> {modes[1]}")
        if name == "candidate" and modes != [["dark"], ["dark"]]:
            ok = False
            print("  FAIL: dark mode state leaks between visitors")
        await a.close()
        await b.close()
    return ok


async def check_perf():
    """Report server time per plot, per gene change, and how long another visitor waits."""
    for name, url in APPS:
        print(f"{name}:")
        for w in WIDGETS:
            c = Client()
            await c.open(url, **{f".clientdata_output_{x}_hidden": x != w for x in WIDGETS})
            times = [(await c.send("update", {"gene": g}))["secs"] for g in ["Gfap", "Snap25", "Aqp4", "Syt1"]]
            print(f"  {w} alone: " + " ".join(f"{t:.2f}" for t in times) + " s")
            await c.close()
        c = Client()
        await c.open(url)
        times = [(await c.send("update", {"gene": g}))["secs"]
                 for g in ["Gfap", "Rbfox3", "P2ry12", "Aldh1l1", "Snap25", "Tmem119"]]
        clicks = [(await c.send("update", {"sex": s}))["secs"] for s in (["Male"], ["Female", "Male"], ["Female"])]
        print("  gene change, all 5 plots visible: " + " ".join(f"{t:.2f}" for t in times)
              + " s | checkbox click: " + " ".join(f"{t:.2f}" for t in clicks) + " s")
        await c.close()
        # A second visitor's trivial update while the first visitor changes gene
        a, b = Client(), Client()
        await a.open(url)
        await b.open(url, gene="Gfap")
        task = asyncio.create_task(a.send("update", {"gene": "Snap25"}))
        await asyncio.sleep(0.05)
        rb = await b.send("update", {"sex": ["Male"]})
        await task
        print(f"  second visitor's update waited {rb['secs']:.2f} s while another visitor changed gene")
        await a.close()
        await b.close()
    return True


CHECKS = {"downloads": check_downloads, "inputs": check_inputs, "darkmode": check_darkmode, "perf": check_perf}

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in CHECKS:
        sys.exit(f"usage: {sys.argv[0]} {'|'.join(CHECKS)}")
    sys.exit(0 if asyncio.run(CHECKS[sys.argv[1]]()) else 1)
