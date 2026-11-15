#!/bin/bash
# Minimum weight of R-round differential characteristics of FeatherCrypt (r = 17) with STP.
# For each R, the model with weight <= W(R) must have a solution and the model with weight <= W(R)-1 must not.
# Usage: ./run_stp.sh [first_round] [last_round]      (default: 2 15)
# STP prints "Valid." when the model has no solution, and a list of ASSERT lines (a characteristic) otherwise.
cd "$(dirname "$0")" || exit 1
STP=${STP:-stp}
PY=${PYTHON:-python3}
W=(0 2 4 8 12 16 19 23 27 31 34 38 42 46 49 53)
mkdir -p work
for R in $(seq "${1:-2}" "${2:-15}"); do
  w=${W[$R]}
  for b in $((w - 1)) "$w"; do
    f=work/FeatherCrypt_${R}r_w${b}.cvc
    $PY FeatherCrypt_stp.py -r "$R" -w "$b" -o "$f" > /dev/null
    if $STP "$f" | head -1 | grep -q "^Valid"; then res="no characteristic"; else res="characteristic found"; fi
    echo "R=$R weight<=$b: $res"
  done
done
