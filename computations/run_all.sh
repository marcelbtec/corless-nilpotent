#!/bin/sh
# Reruns every numerical check. Needs python3 (standard library only; make_figure.py also needs
# matplotlib) and a C compiler with 128-bit integers (gcc or clang).
# The exhaustive search at m = 63 runs 14 slices in parallel; set SLICES to change that.
set -e
cd "$(dirname "$0")"
SLICES=${SLICES:-14}
# run a check and stop if it fails or prints False
run() { out=$("$@") || exit 1; printf '%s\n' "$out"; if printf '%s\n' "$out" | grep -q False; then echo "FAILED: $*"; exit 1; fi; }
run python3 brute.py
run python3 check_proof.py
run python3 check_matrices.py
run python3 mod4.py
run python3 halves.py
run python3 sfraction.py
run python3 check_extras.py
cc -O3 -o enum enum.c
for K in 3 4 5 6; do
  python3 gen_params.py $K > params_$K.txt
  i=0; : > sol_$K.txt; : > log_$K.txt
  while [ $i -lt $SLICES ]; do ./enum params_$K.txt $i $SLICES >> sol_$K.part.$i 2>> log_$K.part.$i & i=$((i+1)); done
  wait
  i=0; while [ $i -lt $SLICES ]; do cat sol_$K.part.$i >> sol_$K.txt; cat log_$K.part.$i >> log_$K.txt; rm sol_$K.part.$i log_$K.part.$i; i=$((i+1)); done
done
run python3 check_enum.py
