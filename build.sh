#!/usr/bin/env bash
# Build lectures.pdf (all twelve lectures, contents page, continuous page numbers).
set -u
cd "$(dirname "$0")"
for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error lectures.tex >/dev/null || { echo "FAIL"; grep -m3 -A3 '^!' lectures.log; exit 1; }; done
echo "ok   lectures.pdf ($(pdfinfo lectures.pdf | awk '/Pages/{print $2}') pages)"
grep -i 'warning' lectures.log | grep -v 'Overfull\|infwarerr' | sed 's/^/     /'
rm -f lectures.aux lectures.out lectures.toc
