#!/bin/bash
# Minimum weight of R-round differential characteristics of FeatherCrypt (r = 17) with STP.
# For each R, the model "weight <= W(R)" must have a solution and the model "weight <= W(R)-1" must not.
# Usage: ./run_stp.sh [first_round] [last_round]      (default: 2 16)
# STP ends its output with "Valid." when the model has no solution and with "Invalid." when it has one
# (the solution, a characteristic, is printed before as ASSERT lines).
# Any other outcome (solver missing, error, time-out, unknown) stops the script with an error.
cd "$(dirname "$0")" || exit 1
STP=${STP:-stp}
PY=${PYTHON:-python3}
command -v "$STP" > /dev/null || { echo "error: STP not found (set STP=/path/to/stp)" >&2; exit 1; }
W=(0 2 4 8 12 16 19 23 27 31 34 38 42 46 49 53 57)
mkdir -p work
"$STP" --version 2>/dev/null | head -1
for R in $(seq "${1:-2}" "${2:-16}"); do
  w=${W[$R]}
  for b in $((w - 1)) "$w"; do
    f=work/FeatherCrypt_${R}r_w${b}.cvc
    "$PY" FeatherCrypt_stp.py -r "$R" -w "$b" -o "$f" > /dev/null || { echo "error: generator failed for $f" >&2; exit 1; }
    "$STP" "$f" > "$f.out" 2> "$f.err"
    rc=$?
    last=$(tail -n 1 "$f.out")
    if [ $rc -ne 0 ]; then echo "error: STP exit code $rc on $f" >&2; exit 1; fi
    case "$last" in
      Valid.*)   res="no characteristic" ;;
      Invalid.*) res="characteristic found" ;;
      *) echo "error: unexpected STP output on $f: '$last'" >&2; exit 1 ;;
    esac
    echo "R=$R weight<=$b: $res"
    if [ "$b" -eq "$w" ] && [ "$res" != "characteristic found" ]; then echo "error: expected a characteristic" >&2; exit 1; fi
    if [ "$b" -lt "$w" ] && [ "$res" != "no characteristic" ]; then echo "error: expected no characteristic" >&2; exit 1; fi
  done
done
echo "all checks passed"
