#!/bin/bash
# Check the candidate app's page: plotly.js and the output binding load (in order, after
# Shiny), are served correctly, and nothing is loaded from external hosts.
# Run via tests/run.sh (uses CAND_URL and WORK).
CAND_URL=${CAND_URL:-http://127.0.0.1:8102}
WORK=${WORK:-${TMPDIR:-/tmp}/aging-app-tests}
REPO=$(cd "$(dirname "$0")/.." && pwd)
page=$WORK/candidate_page.html
status=0

curl -s "$CAND_URL/" > "$page"
echo "script tags, in order:"
grep -o '<script[^>]*src="[^"]*"' "$page" | sed 's/.*src="//; s/"$//' | nl

shiny=$(grep -o '<script[^>]*src="[^"]*"' "$page" | grep -n 'shiny.js"' | cut -d: -f1)
plotly=$(grep -o '<script[^>]*src="[^"]*"' "$page" | grep -n 'plotly.min.js"' | cut -d: -f1)
binding=$(grep -o '<script[^>]*src="[^"]*"' "$page" | grep -n 'plotly_output.js"' | cut -d: -f1)
if [[ -z $shiny || -z $plotly || -z $binding || ! ($shiny -lt $plotly && $plotly -lt $binding) ]]; then
  echo "FAIL: expected shiny.js, then plotly.min.js, then plotly_output.js"; status=1
fi

for asset in $(grep -oE '(src|href)="[^"]*plotly[^"]*"' "$page" | sed 's/.*="//; s/"$//'); do
  code=$(curl -s -o "$WORK/asset" -w '%{http_code}' "$CAND_URL/$asset")
  echo "GET $asset -> $code, $(wc -c < "$WORK/asset") bytes"
  [[ $code == 200 ]] || { echo "FAIL: $asset not served"; status=1; }
  case $asset in
    */plotly_output.js|*/plotly_output.css)
      cmp -s "$WORK/asset" "$REPO/www/$(basename "$asset")" || { echo "FAIL: $asset differs from www/"; status=1; };;
  esac
done

external=$(grep -oE '(src|href|data-[a-z-]*cdn)="https?://[^"]*"' "$page" | sort -u)
if [[ -n $external ]]; then
  echo "FAIL: page references external URLs: $external"; status=1
else
  echo "no external URLs referenced"
fi
exit $status
