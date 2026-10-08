#!/bin/bash
# Start a reference version of the app (from git) and the working tree ("candidate") side by
# side on localhost, run the requested checks against both, then stop them.
#
# usage: tests/run.sh [check ...]
#   checks: plots downloads inputs darkmode page perf screenshots
#   default: plots downloads inputs darkmode page perf   (screenshots needs Chrome)
#
# environment:
#   REF_REV      git revision used as the reference              (default: HEAD)
#   CAND_PYTHON  python of a venv with requirements.txt          (default: python3)
#   REF_PYTHON   python of a venv with the reference's requirements (default: CAND_PYTHON)
#   TEST_PYTHON  python for the checks; screenshots also need tests/requirements.txt
#                                                                 (default: CAND_PYTHON)
#   REF_PORT, CAND_PORT  localhost ports                          (default: 8101, 8102)
#   WORK         scratch directory for the reference copy, logs and screenshots
#                                                                 (default: $TMPDIR/aging-app-tests)
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd)
REF_REV=${REF_REV:-HEAD}
CAND_PYTHON=${CAND_PYTHON:-python3}
REF_PYTHON=${REF_PYTHON:-$CAND_PYTHON}
TEST_PYTHON=${TEST_PYTHON:-$CAND_PYTHON}
REF_PORT=${REF_PORT:-8101}
CAND_PORT=${CAND_PORT:-8102}
export WORK=${WORK:-${TMPDIR:-/tmp}/aging-app-tests}
export REF_URL=http://127.0.0.1:$REF_PORT CAND_URL=http://127.0.0.1:$CAND_PORT
export PYTHONDONTWRITEBYTECODE=1  # keep __pycache__ out of the repo
checks=("$@")
[[ ${#checks[@]} -eq 0 ]] && checks=(plots downloads inputs darkmode page perf)

# Reference: the app code at REF_REV, using the working tree's data files
ref_dir=$WORK/reference
rm -rf "$ref_dir" && mkdir -p "$ref_dir"
paths=$(git -C "$REPO" ls-tree --name-only "$REF_REV" app.py www) || exit 1
git -C "$REPO" archive "$REF_REV" $paths | tar -x -C "$ref_dir" || exit 1
for f in "$REPO"/*.parquet "$REPO"/DETECTED_GENES.csv; do ln -s "$f" "$ref_dir/"; done
echo "reference: $(git -C "$REPO" log --oneline -1 "$REF_REV")"
echo "candidate: working tree of $REPO"

(cd "$ref_dir" && exec "$REF_PYTHON" -m shiny run --host 127.0.0.1 --port "$REF_PORT" app.py) > "$WORK/reference.log" 2>&1 &
ref_pid=$!
(cd "$REPO" && exec "$CAND_PYTHON" -m shiny run --host 127.0.0.1 --port "$CAND_PORT" app.py) > "$WORK/candidate.log" 2>&1 &
cand_pid=$!
trap 'kill $ref_pid $cand_pid 2>/dev/null; wait 2>/dev/null' EXIT

for _ in $(seq 120); do
  curl -sf -o /dev/null "$REF_URL" && curl -sf -o /dev/null "$CAND_URL" && break
  kill -0 $ref_pid 2>/dev/null && kill -0 $cand_pid 2>/dev/null || break
  sleep 1
done
if ! curl -sf -o /dev/null "$REF_URL" || ! curl -sf -o /dev/null "$CAND_URL"; then
  echo "apps did not start; see $WORK/reference.log and $WORK/candidate.log"; exit 1
fi
echo "memory at startup: reference $(ps -o rss= -p $ref_pid) KB, candidate $(ps -o rss= -p $cand_pid) KB"

cd "$REPO"
failed=()
for check in "${checks[@]}"; do
  echo "== $check"
  case $check in
    plots) "$TEST_PYTHON" tests/compare_plots.py;;
    downloads|inputs|darkmode|perf) "$TEST_PYTHON" tests/harness.py "$check";;
    page) bash tests/pagecheck.sh;;
    screenshots) "$TEST_PYTHON" tests/screenshots.py;;
    *) echo "unknown check: $check"; false;;
  esac || failed+=("$check")
done
echo "memory after checks: reference $(ps -o rss= -p $ref_pid) KB, candidate $(ps -o rss= -p $cand_pid) KB"

if [[ ${#failed[@]} -gt 0 ]]; then
  echo "FAILED: ${failed[*]}  (server logs in $WORK)"; exit 1
fi
echo "all checks passed: ${checks[*]}"
