#!/usr/bin/env bash
# Build latexdiff PDFs comparing the current files in the repo root with the
# proposed revisions in revised/.  Output: diffs/<name>_diff.pdf (one per file)
# and diffs/all_diffs.pdf (everything concatenated, in course order).
#
# Usage:  ./make_diffs.sh            # all files
#         ./make_diffs.sh lecture3   # one file
set -u
cd "$(dirname "$0")"
mkdir -p diffs
ln -sfn ../revised/figures diffs/figures
ln -sfn ../mycommands.sty diffs/mycommands.sty

if [ $# -gt 0 ]; then
  files=("$@")
else
  files=(lecture1 lecture2 lecture3 lecture4 lecture5 lecture6 lecture7 lecture8 \
         lecture9 lecture10 lecture11 lecture12 problem1 solution1 problem2 solution2)
fi

ok=()
for f in "${files[@]}"; do
  if [ ! -f "revised/$f.tex" ]; then echo "skip $f (no revised file)"; continue; fi
  if [[ "$f" == lecture* ]]; then
    # revised lectures are bodies only: wrap in the original's preamble so latexdiff sees two full documents
    sed -n '1,/\\maketitle/p' "$f.tex" > "diffs/${f}_new.tex"
    printf '\\input{../revised/%s}\n\\end{document}\n' "$f" >> "diffs/${f}_new.tex"
    latexdiff --flatten --type=UNDERLINE --math-markup=whole "$f.tex" "diffs/${f}_new.tex" \
        > "diffs/${f}_diff.tex" 2> "diffs/${f}_diff.latexdiff.log"
    rm -f "diffs/${f}_new.tex"
  else
    latexdiff --type=UNDERLINE --math-markup=whole "$f.tex" "revised/$f.tex" \
        > "diffs/${f}_diff.tex" 2> "diffs/${f}_diff.latexdiff.log"
  fi
  ( cd diffs && pdflatex -interaction=nonstopmode "${f}_diff.tex" > /dev/null 2>&1 \
             && pdflatex -interaction=nonstopmode "${f}_diff.tex" > /dev/null 2>&1 )
  if [ -f "diffs/${f}_diff.pdf" ] && ! grep -q '^!' "diffs/${f}_diff.log"; then
    echo "ok   $f  ($(grep -c 'DIFadd\|DIFdel' "diffs/${f}_diff.tex") markup commands, $(pdfinfo "diffs/${f}_diff.pdf" | awk '/Pages/{print $2}') pages)"
    ok+=("diffs/${f}_diff.pdf")
  else
    echo "FAIL $f  -- see diffs/${f}_diff.log"
    grep -m3 '^!' "diffs/${f}_diff.log" || true
  fi
done

if [ ${#ok[@]} -gt 1 ] && command -v pdfunite > /dev/null; then
  pdfunite "${ok[@]}" diffs/all_diffs.pdf && echo "wrote diffs/all_diffs.pdf"
fi
