#!/bin/bash
# Clone a generated run directory once per value of a GIZMO parameter, for a parameter survey.
# Usage: make_grid.sh <run dir> <parameter> <value> [<value> ...]
#   e.g. make_grid.sh runs/mw SfEffPerFreeFall 1 0.1 0.03 0.01   ->  runs/mw-SfEffPerFreeFall1, ...
# The run directory must hold the binary as GIZMO. ICs and binary are hard-linked where the filesystem allows.
set -e
src=${1%/}; par=$2; shift 2
grep -q "^$par " "$src/mw.params" || { echo "$par is not in $src/mw.params"; exit 1; }
for v in "$@"; do
    dst=$src-$par$v
    mkdir "$dst"
    for f in mw_ics.hdf5 GIZMO; do ln "$src/$f" "$dst/$f" 2> /dev/null || cp "$src/$f" "$dst/$f"; done
    cp "$src/Config.sh" "$dst/"
    for f in TREECOOL spcool_tables; do [ -e "$src/$f" ] && cp -r "$src/$f" "$dst/"; done
    sed -E "s|^($par +)[^ ]+|\1$v|" "$src/mw.params" > "$dst/mw.params"
    echo "$dst: $(grep "^$par " "$dst/mw.params" | tr -s ' ')"
done
