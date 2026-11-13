#!/usr/bin/env bash
# Build demo/DEMO.html from demo/DEMO.md, for reading the demo offline in a
# browser after a clone or ZIP download (GitHub renders the .md itself).
#
#   bash scripts/build_demo.sh
#
# Same pandoc setup as scripts/build_userguide.sh, with two differences:
#   - figures are NOT embedded: the HTML loads them from demo/img/, which ships
#     beside it, so image changes do not re-add megabytes to git history;
#   - links to the User Guide .md are rewritten to its .html copy.
# The stylesheet is inlined, so the page needs nothing outside demo/ and docs/.
set -euo pipefail
cd "$(dirname "$0")/.."

SRC="demo/DEMO.md"
OUT="demo/DEMO.html"
CSS="scripts/userguide.css"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# .md links that a browser cannot render -> their HTML / web equivalents
sed -e 's#(\.\./docs/UserGuide\.v1\.0\.md#(../docs/UserGuide.v1.0.html#g' \
    -e 's#(\.\./README\.md)#(https://github.com/weinstein-bioinfo/EMtrace\#readme)#g' \
    "$SRC" > "$tmp/DEMO.md"

# Light mode only: drop the stylesheet's dark-mode block (its light text would
# otherwise sit on the page background) and declare the page light.
{ echo "<style>"
  echo ":root { color-scheme: light; }"
  awk '/@media \(prefers-color-scheme: dark\)/ {skip=1} !skip {print} skip && /^}/ {skip=0}' "$CSS"
  echo "</style>"; } > "$tmp/header.html"

# Run from demo/ so the img/ paths resolve exactly as they do on GitHub
# (gfm without auto heading ids: DEMO.md sets an explicit <a id> per section,
# and pandoc's own ids would duplicate them)
(cd demo && pandoc "$tmp/DEMO.md" -f gfm-gfm_auto_identifiers -t html5 --standalone \
    --metadata pagetitle="EMtrace tutorial" \
    -M document-css=false \
    --include-in-header "$tmp/header.html" \
    -o "DEMO.html")

echo "Wrote $OUT"
