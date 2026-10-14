# App comparison tests

These tests start two versions of the app side by side on localhost and drive them the way a
browser does, over Shiny's websocket protocol:

- **reference**: the app code at a git revision (`REF_REV`, default `HEAD`), using the
  working tree's data files
- **candidate**: the working tree

They answer "did my change alter what users see, or how fast it is?"

| Check | What it does | Fails when |
|---|---|---|
| `plots` | Captures what each app sends the browser for every plot (trace data, layout including the template, plotly config) over 4 genes, all 7 groupings, filter subsets, and light and dark mode. Floats are compared exactly. | any plot differs |
| `downloads` | Compares the three CSV downloads over 3 genes × 7 groupings × 2 filter sets | any download differs |
| `inputs` | Sends forged inputs (unknown gene, wrong types, invalid grouping) and a gene with no methylation data | the candidate shows raw exception text |
| `darkmode` | Two visitors each click the dark-mode toggle once | one visitor's toggle affects the other's |
| `genesearch` | Searches genes through the Gene box's server-side search (e.g. "gfa" finds Gfap) | the search doesn't work or the gene list is written into the page |
| `page` | Checks the page HTML size and that plotly.js and `www/plotly_output.js`/`.css` load after Shiny and are served correctly | page over 200 KB, assets missing or out of order, or any external URL on the page |
| `perf` | Reports server time per plot, per gene change, and how long a second visitor waits | (report only) |
| `screenshots` | Renders both apps in headless Chrome, measures page load (time to first plot, main-thread blocking), picks genes by typing over the current gene in the Gene box, checks that leaving the box without picking keeps the gene, and compares every plot pixel by pixel in 5 scenarios, saving images to `$WORK/shots/` | the candidate's page load blocks the browser for over 2 s at once, typing over the gene doesn't work, leaving the box changes the gene, or it has a JavaScript error, a missing or squashed plot, or a correlation plot drawn narrower than its container |

`plots` handles both ways the app has sent figures (shinywidgets widgets and plain JSON), so the
reference can be a revision from before the switch to plotly.js.

## Running

Make a venv with the app's requirements, plus `tests/requirements.txt` if you'll run
`screenshots`:

`requirements.txt` is locked for Python 3.14 (the container's Python), with hashes, so the
venv needs Python 3.14 (on the OMRF cluster: `module load python/3.14.7`). Install the
test-only packages in a separate step, because pip requires hashes for everything installed
together with a hash-locked file.

The cluster's python modules set `PYTHONPATH` to the module's own site-packages, which comes
ahead of the venv: pip skips packages it finds there, and the app imports the module's numpy,
scipy and pandas instead of the locked ones. Unset it before installing or running the app
(`run.sh` unsets it itself):

```bash
unset PYTHONPATH
python3.14 -m venv /s/ai/$USER/venvs/aging-app
/s/ai/$USER/venvs/aging-app/bin/pip install --require-hashes -r requirements.txt
/s/ai/$USER/venvs/aging-app/bin/pip install -r tests/requirements.txt   # for screenshots

CAND_PYTHON=/s/ai/$USER/venvs/aging-app/bin/python tests/run.sh               # default checks
CAND_PYTHON=/s/ai/$USER/venvs/aging-app/bin/python tests/run.sh screenshots   # needs Chrome
```

To test an app that is already running, such as the container, set `CAND_URL` instead of
starting the working tree:

```bash
CAND_URL=http://127.0.0.1:8000 CAND_PYTHON=/s/ai/$USER/venvs/aging-app/bin/python tests/run.sh
```

If the reference revision has different requirements, give it its own venv with `REF_PYTHON`
(for example, for revisions that still use shinywidgets, install that revision's
`requirements.txt`: `git show <rev>:requirements.txt > ref-requirements.txt`).

Other settings: `REF_REV`, `TEST_PYTHON`, `REF_PORT`/`CAND_PORT` (default 8101/8102),
`WORK` (scratch directory, default `$TMPDIR/aging-app-tests`; holds the reference copy,
server logs, and screenshots), `CHROME` (default `/usr/bin/google-chrome`).

The two apps take about 0.6 GB of memory each, and the default checks take a few minutes. On
the cluster, run them in an allocation, for example:

```bash
srun --cpus-per-task=2 --mem=4G --time=00:30:00 tests/run.sh
```

`screenshots` needs Chrome on the machine running it, and Chrome won't start inside Claude
Code's sandbox (it isn't allowed to create the sockets it needs).

## Baseline results

Reference `259c6d6` (security fixes, plots still drawn with shinywidgets) vs. `4645b3a` (plots
drawn with plotly.js directly), 2026-10-08:

| | shinywidgets (`259c6d6`) | plotly.js (`4645b3a`) |
|---|---|---|
| `plots`: 165 renders compared | | identical |
| `downloads`: 126 compared | | identical |
| Server time, gene change with all 5 plots | ~3.3 s | ~0.3 s |
| Server time per plot | 0.5–0.75 s | 0.01–0.09 s |
| A second visitor's wait while another changes gene | ~7 s | ~0.4 s |
| Data sent per plot render | 5.3 MB | 11 KB |
| Memory after the checks | ~790 MB | ~570 MB |

`screenshots`: 19 of 25 plot screenshots were pixel-identical. The other 6 are correlation plots
that shinywidgets drew at plotly's 700px default width after being redrawn while their tab was
hidden; `4645b3a` draws them full width. The "No data" messages for a gene without methylation
data are 48px tall in both versions (expected; the width check flags them as narrow but doesn't
fail on them).

Gene box, `54a78e5` (all 25,114 genes written into the page) vs. server-side gene search,
2026-10-09, headless Chrome: page load to first plot ~27 s → ~2.2 s; longest main-thread block
~25 s → ~1 s; page HTML 1.5 MB → 13 KB. Typing a gene shows the suggestion in about 1 s in both.
All 25 plot screenshots were pixel-identical.
