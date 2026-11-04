#!/usr/bin/env bash
# Runs SOLARIS on the GaMD trajectory of each enzyme, one at a time, each in its own folder.
#   bash tools/run_gamd_per_enzyme.sh               (all enzymes)
#   bash tools/run_gamd_per_enzyme.sh ScTYR BmTYR   (only these)
set -u
CONFIG=examples/TYR_GaMD.json
ENZYMES=("$@")
[ ${#ENZYMES[@]} -eq 0 ] && ENZYMES=(ScTYR VvTYR AbTYR BmTYR HsTYR VsTYR JrTYR)
for ENZ in "${ENZYMES[@]}"; do
    echo "=================== $ENZ ==================="
    python3 solaris.py "$CONFIG" --only "$ENZ" --states GaMD --out "SOLARIS_GaMD/$ENZ" 2>&1 | tee "SOLARIS_GaMD_${ENZ}.log"
done
