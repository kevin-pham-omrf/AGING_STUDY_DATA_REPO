"""Render the reference and candidate apps in headless Chrome and compare screenshots of
every plot, pixel by pixel, across several scenarios (genes, groupings, dark mode, a gene
without methylation data), picking genes by typing in the Gene box. Also records page
load time and main-thread blocking, layout details for missing or squashed plots, browser
console errors, and the drawn width of the correlation plots.

Needs Chrome (CHROME, default /usr/bin/google-chrome) and tests/requirements.txt.
Results go to $WORK/shots/ (report.txt plus PNGs; DIFF__*.png are side-by-side
reference | candidate | changed pixels). Run via tests/run.sh.

Fails if the candidate's page load blocks the browser for over 2 s at once, or it has a
JavaScript error, a plot that is missing or squashed (other
than the expected "no data" messages), or a correlation plot drawn narrower than its
container. Pixel differences from the reference are reported, not failed: some are
expected (see tests/README.md).
"""
import os
import sys
import time
from pathlib import Path

from PIL import Image, ImageChops
from playwright.sync_api import TimeoutError as PWTimeout
from playwright.sync_api import sync_playwright

APPS = {"reference": os.environ.get("REF_URL", "http://127.0.0.1:8101"),
        "candidate": os.environ.get("CAND_URL", "http://127.0.0.1:8102")}
TABS = {
    "Gene Expression": ["expression_plot"],
    "Gene DNA Modifications": ["tss_plot", "gene_body_plot"],
    "Correlation Plots": ["tss_corr_plot", "gene_corr_plot"],
}
NO_DATA_SCENARIO = "no_methylation_gene_light"
# (name, gene to pick by typing in the Gene box, JS run before the screenshots) --
# applied in order, state carries over
SCENARIOS = [
    ("default", None, None),
    ("Gfap_age+sex", "Gfap", "$('#filter').val('4').trigger('change');"),
    ("Gfap_sex_dark", None, "$('#filter').val('2').trigger('change'); $('#toggle_dark').click();"),
    ("Snap25_celltype_subset_dark", "Snap25", "$('#filter').val('3').trigger('change');"
                                              "$('input[name=line][value=Microglia]').click();"),
    (NO_DATA_SCENARIO, "4732440D04Rik", "$('#toggle_dark').click(); $('#filter').val('7').trigger('change');"),
]
WORK = Path(os.environ.get("WORK", Path(os.environ.get("TMPDIR", "/tmp")) / "aging-app-tests"))
OUT = WORK / "shots"
OUT.mkdir(parents=True, exist_ok=True)
for old in OUT.glob("*.png"):
    old.unlink()
report = open(OUT / "report.txt", "w")
failures = []


def log(msg):
    print(msg)
    report.write(msg + "\n")
    report.flush()


def fail(msg):
    failures.append(msg)
    log(f"FAIL: {msg}")


# Element and its ancestors: size and the CSS that drives fill layout
BOX_CHAIN = """(id) => { const out = []; let el = document.getElementById(id);
  for (let i = 0; el && i < 14; i++, el = el.parentElement) {
    const cs = getComputedStyle(el), r = el.getBoundingClientRect();
    out.push(`${el.tagName.toLowerCase()}${el.id ? '#' + el.id : ''}.${[...el.classList].join('.')}`
      + ` h=${Math.round(r.height)} top=${Math.round(r.top)} display=${cs.display} flex=${cs.flex}`
      + ` minH=${cs.minHeight} overflow=${cs.overflow} vis=${cs.visibility}`);
  } return out; }"""


def pick_gene(page, app, gene):
    """Type into the Gene box and click the matching suggestion, as a visitor would.

    Logs how long the suggestion took to appear. If it never does, saves a screenshot and
    selects the gene directly so the remaining scenarios can still run.
    """
    start = time.perf_counter()
    page.locator("#gene + .selectize-control .selectize-input").click()
    # While a gene is selected the box is "full" and selectize ignores typed characters;
    # Backspace removes the current gene first. Key presses go to selectize's text field.
    field = page.locator("#gene + .selectize-control .selectize-input input")
    field.press("Backspace")
    field.press_sequentially(gene, delay=50)
    try:
        page.locator(f"#gene + .selectize-control .selectize-dropdown .option[data-value='{gene}']").click(timeout=120000)
        log(f"[{app}] typed {gene!r} in the Gene box: suggestion clicked after {time.perf_counter() - start:.1f} s")
    except PWTimeout:
        page.screenshot(path=OUT / f"PICK_FAILED__{app}__{gene}.png", full_page=True)
        msg = f"typed {gene!r} in the Gene box: no suggestion after {time.perf_counter() - start:.0f} s"
        if app == "candidate":
            fail(msg)
        else:
            log(f"[{app}] {msg}; selecting it directly instead")
        page.keyboard.press("Escape")
        page.evaluate("(g) => { const s = $('#gene')[0].selectize; s.addOption({value: g, label: g}); s.setValue(g); }", gene)


# Main-thread tasks over 50 ms (what makes a page unresponsive), recorded from page start
LONG_TASKS = """window.__longTasks = [];
new PerformanceObserver(l => l.getEntries().forEach(e => window.__longTasks.push(e.duration)))
  .observe({type: 'longtask', buffered: true});"""


def settle(page):
    page.wait_for_timeout(400)
    page.wait_for_function("() => !document.documentElement.classList.contains('shiny-busy')", timeout=60000)
    page.wait_for_timeout(3000)  # client-side plot rendering


def on_page_error(app, msg):
    log(f"[{app}] {msg}")
    if app == "candidate":
        failures.append(msg)


with sync_playwright() as p:
    # Chrome writes crash reports etc. under $HOME; keep that in the work directory
    home = str(WORK / "chromehome")
    browser = p.chromium.launch(executable_path=os.environ.get("CHROME", "/usr/bin/google-chrome"), headless=True,
                                args=["--no-sandbox", "--disable-crash-reporter", "--disable-gpu"],
                                env={**os.environ, "HOME": home, "XDG_CONFIG_HOME": home, "XDG_CACHE_HOME": home})
    for app, url in APPS.items():
        page = browser.new_page(viewport={"width": 1500, "height": 1100})
        page.add_init_script(LONG_TASKS)
        page.on("pageerror", lambda e, app=app: on_page_error(app, f"PAGE ERROR: {e}"))
        page.on("console", lambda m, app=app: m.type == "error" and log(f"[{app}] console error: {m.text[:300]}"))
        start = time.perf_counter()
        page.goto(url, timeout=600000)
        page.wait_for_selector("#expression_plot .main-svg", timeout=600000)
        load_secs = time.perf_counter() - start
        tasks = page.evaluate("() => window.__longTasks")
        log(f"[{app}] page load to first plot: {load_secs:.1f} s; main thread blocked {sum(tasks) / 1000:.1f} s "
            f"in {len(tasks)} long tasks, longest {max(tasks, default=0) / 1000:.1f} s")
        if app == "candidate" and max(tasks, default=0) > 2000:
            fail(f"page load blocked the browser for {max(tasks) / 1000:.1f} s at once")
        settle(page)
        for name, gene, js in SCENARIOS:
            if gene:
                pick_gene(page, app, gene)
            if js:
                page.evaluate(f"() => {{ {js} }}")
            for tab, plots in TABS.items():
                page.get_by_role("tab", name=tab).click()
                settle(page)
                page.screenshot(path=OUT / f"PAGE__{app}__{name}__{tab.replace(' ', '_')}.png", full_page=True)
                for plot in plots:
                    box = page.locator(f"#{plot}").bounding_box()
                    if not box or box["height"] < 100:
                        log(f"[{app}] {name} / {tab} / {plot}: missing or squashed, box={box}")
                        for line in page.evaluate(BOX_CHAIN, plot):
                            log(f"      {line}")
                        if app == "candidate" and not (name == NO_DATA_SCENARIO and tab != "Gene Expression"):
                            failures.append(f"{name} / {plot} missing or squashed")
                    try:
                        page.locator(f"#{plot}").screenshot(path=OUT / f"{app}__{name}__{plot}.png", timeout=5000)
                    except PWTimeout:
                        log(f"[{app}] {name} / {plot}: no screenshot (element not visible)")
        page.close()
    browser.close()

identical = total = 0
for f in sorted(OUT.glob("reference__*.png")):
    label = f.stem.replace("reference__", "")
    other = OUT / f.name.replace("reference__", "candidate__")
    total += 1
    if not other.exists():
        log(f"  {label}: no candidate screenshot")
        continue
    a, b = Image.open(f).convert("RGB"), Image.open(other).convert("RGB")
    if a.size != b.size:
        log(f"  {label}: SIZE differs {a.size} vs {b.size}")
        continue
    d = ImageChops.difference(a, b)
    n = sum(1 for px in d.getdata() if px != (0, 0, 0))
    if n == 0:
        identical += 1
    else:
        log(f"  {label}: {n} of {a.size[0] * a.size[1]} pixels differ, region {d.getbbox()}")
        side = Image.new("RGB", (a.size[0] * 3, a.size[1]), "white")
        side.paste(a, (0, 0))
        side.paste(b, (a.size[0], 0))
        side.paste(d.point(lambda v: 255 if v else 0), (a.size[0] * 2, 0))
        side.save(OUT / f"DIFF__{label}.png")
log(f"screenshots compared: {total}, pixel-identical: {identical}  (images in {OUT})")


def right_edge(path):
    """x of the rightmost drawn pixel (compared with the background color)."""
    im = Image.open(path).convert("RGB")
    w, h = im.size
    px = im.load()
    bg = px[w - 1, h // 2]
    return max((x for x in range(w) for y in range(0, h, 4) if px[x, y] != bg), default=0), w


log("correlation plot widths (where drawn content ends; full width is ~98% of the image):")
for f in sorted(OUT.glob("*corr_plot.png")):
    if f.name.startswith(("reference__", "candidate__")):
        edge, w = right_edge(f)
        narrow = edge < 0.8 * w
        log(f"  {f.stem}: {edge} of {w}{'  <-- NARROW' if narrow else ''}")
        # The "no data" scenario shows a short text message instead of a plot
        if narrow and f.name.startswith("candidate__") and NO_DATA_SCENARIO not in f.name:
            failures.append(f"{f.stem} drawn narrow")

log(f"candidate failures: {len(failures)}" + "".join(f"\n  {m}" for m in failures))
sys.exit(1 if failures else 0)
