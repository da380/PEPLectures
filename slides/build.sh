#!/usr/bin/env bash
# Build the lecture slides: ./build.sh [13 14 ...]  (default: every slidesNN.tex present).
# slidesNN.pdf is the plain version, for uploading. A deck with animations (\animation) is also built as
# slidesNN-live.pdf, in which their stills open the movies in figures/anim when clicked; present from slides/.
set -u
cd "$(dirname "$0")"
if [ $# -gt 0 ]; then nums=("$@"); else nums=($(ls slides[0-9][0-9].tex | sed 's/slides\([0-9]*\)\.tex/\1/')); fi
run() {   # run <jobname> <tex input>
  local job=$1 input=$2 ok=1
  for i in 1 2; do
    pdflatex -interaction=nonstopmode -halt-on-error -jobname "$job" "$input" >/dev/null || { ok=0; break; }
  done
  if [ $ok = 0 ]; then echo "FAIL $job"; grep -m3 -A3 '^!' $job.log; return; fi
  echo "ok   $job.pdf ($(pdfinfo $job.pdf | awk '/Pages/{print $2}') slides)"
  grep -E 'Overfull|Warning' $job.log | grep -v 'infwarerr' | sed 's/^/     /'
  rm -f $job.aux $job.nav $job.out $job.snm $job.toc $job.vrb
}
for n in "${nums[@]}"; do
  f=slides$n
  run $f $f.tex
  if grep -q '\\animation' $f.tex; then run $f-live "\def\peplive{}\input{$f.tex}"; fi
done
