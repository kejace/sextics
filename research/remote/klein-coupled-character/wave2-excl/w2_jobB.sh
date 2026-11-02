#!/bin/bash
# Job B (remote): msolve Groebner runs (leading ideal, grevlex, GF(p)) for the full (non-eigen) Hessian-relaxation systems.
# usage: bash w2_jobB.sh <threads> <timeout-seconds> file1.ms file2.ms ...
T=$1; TO=$2; shift 2
for f in "$@"; do
  echo "$(date +%T) START $f"
  /usr/bin/time -v timeout $TO msolve -f "$f" -o "${f%.ms}.lead" -g 1 -t $T -v 1 > "${f%.ms}.msolve.log" 2>&1
  echo "$(date +%T) END $f rc=$?"
  grep -E "Maximum resident|Elapsed" "${f%.ms}.msolve.log"
done
echo JOBB_DONE
