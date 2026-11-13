#!/usr/bin/env bash
# Render the User Guide to a single self-contained HTML file.
#
# Usage:
#   bash scripts/build_userguide.sh
#
# Output:
#   docs/UserGuide.v1.0.html  (opened by Help > User Guide; the .md stays the
#   source). Commit it, so builds on machines without pandoc ship it too.
#
# Prerequisites:
#   - pandoc 3.x (brew install pandoc)
#
# gfm input gives GitHub-style heading ids, which is what the guide's TOC
# links were written against.

set -euo pipefail
cd "$(dirname "$0")/.."

SRC="docs/UserGuide.v1.0.md"
OUT="docs/UserGuide.v1.0.html"

# Links to other .md pages point to their HTML copies here (a browser shows .md
# as plain text); the .md source keeps .md links for GitHub.
sed -e 's#(\.\./demo/DEMO\.md#(../demo/DEMO.html#g' "$SRC" |
pandoc -f gfm -t html5 --standalone --embed-resources \
    --resource-path=docs \
    --metadata pagetitle="EMtrace User Guide" \
    --css scripts/userguide.css \
    -o "$OUT"

echo "Wrote $OUT"
