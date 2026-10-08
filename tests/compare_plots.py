"""Compare what the reference and candidate apps send the browser to draw each plot:
trace data, layout (template included) and plotly config, with floats compared exactly.

Handles both ways the app has sent figures: shinywidgets widgets (before the switch to
plain plotly.js) and plain JSON figures. Run via tests/run.sh.
"""
import asyncio
import base64
import json
import sys

import numpy as np

from harness import ALL_SELECTED, CAND_URL, REF_URL, WIDGETS, Client, genes_with_all_data


def nan_to_none(values):
    return [None if v != v else v for v in values]


def figures(result):
    """{output: (data, layout, config)} for every plot rendered in this round trip."""
    by_model = {}
    for m in result["custom"]:
        if "shinywidgets_comm_open" not in m:
            continue
        msg = m["shinywidgets_comm_open"]
        msg = json.loads(msg) if isinstance(msg, str) else msg
        state = msg["content"]["data"]["state"]
        # Numeric arrays travel as base64 buffers alongside the JSON state
        for path, buf in zip(msg["content"]["data"]["buffer_paths"], msg["buffers"]):
            node = state
            for key in path[:-2]:
                node = node[key]
            spec = node[path[-2]]
            arr = np.frombuffer(base64.b64decode(buf), dtype=spec["dtype"]).reshape(spec["shape"])
            node[path[-2]] = nan_to_none(arr.tolist())
        for trace in state["_widget_data"]:
            trace.pop("uid", None)  # random per widget
        # plotly.py widget default; only used by the (disabled) "Edit in Chart Studio" link
        state["_config"].pop("plotlyServerURL", None)
        comm_id = msg["content"].get("comm_id") or msg.get("comm_id") or msg["header"]["msg_id"]
        by_model[comm_id] = (state["_widget_data"], state["_widget_layout"], state["_config"])
    out = {}
    for name, v in result["values"].items():
        if isinstance(v, dict) and v.get("model_id") in by_model:
            out[name] = by_model[v["model_id"]]
        elif isinstance(v, dict) and "data" in v:
            out[name] = (v["data"], v["layout"], v.get("config", {}))
    return out


def diff(a, b, path=""):
    """Paths where two JSON-like values differ (floats compared exactly)."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}.{k}: only in {'reference' if k in a else 'candidate'}")
            else:
                out += diff(a[k], b[k], f"{path}.{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in diff(x, y, f"{path}[{i}]")]
    numbers = isinstance(a, (int, float)) and isinstance(b, (int, float))
    return [] if a == b or (numbers and float(a) == float(b)) else [f"{path}: {json.dumps(a)[:60]} vs {json.dumps(b)[:60]}"]


async def main():
    genes = ["Cx3cr1", "Gfap"] + genes_with_all_data(2, seed=7)
    subset = {"age": ["Young", "Old"], "sex": ["Male"], "line": ["Neurons", "Astrocytes"]}
    ref, cand = Client(), Client()
    ref.quiet = cand.quiet = 0.5
    rounds = [("initial", await ref.open(REF_URL), await cand.open(CAND_URL))]
    clicks = 0
    for g in genes:
        for flt in "1234567":
            sel = subset if int(flt) % 2 else {}
            update = {"gene": g, "filter": flt, **ALL_SELECTED, **sel}
            rounds.append((f"{g} grouping {flt}{' subset' if sel else ''}",
                           await ref.send("update", update), await cand.send("update", update)))
        clicks += 1  # flip dark mode once per gene, so both themes get covered
        toggle = {"toggle_dark:shiny.action": clicks}
        rounds.append((f"{g} dark={clicks % 2 == 1}", await ref.send("update", toggle), await cand.send("update", toggle)))

    compared = identical = 0
    sizes = {"reference": [], "candidate": []}
    for label, r_ref, r_cand in rounds:
        f_ref, f_cand = figures(r_ref), figures(r_cand)
        for out in WIDGETS:
            if out not in f_ref and out not in f_cand:
                continue
            if (out in f_ref) != (out in f_cand):
                print(f"  {label} {out}: rendered by only one app ({'reference' if out in f_ref else 'candidate'})")
                compared += 1
                continue
            compared += 1
            d = diff(list(f_ref[out]), list(f_cand[out]))
            if d:
                print(f"  DIFF {label} {out}: {len(d)} differences, e.g. " + "; ".join(d[:4]))
            else:
                identical += 1
        for name, r in (("reference", r_ref), ("candidate", r_cand)):
            sizes[name] += [len(json.dumps(m)) for m in r["custom"] if "shinywidgets_comm_open" in m]
            sizes[name] += [len(json.dumps(v)) for v in r["values"].values() if isinstance(v, dict) and "data" in v]
    print(f"plots compared: {compared}, identical (data + layout + config): {identical}")
    for name, s in sizes.items():
        if s:
            print(f"  bytes sent per plot render, {name}: median {int(np.median(s)):,} (n={len(s)})")
    return compared > 0 and identical == compared


if __name__ == "__main__":
    sys.exit(0 if asyncio.run(main()) else 1)
