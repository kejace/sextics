#!/bin/bash
# usage: run_msolve.sh <input.ms> <threads> <tag>   -- runs msolve (leading ideal + quotient dimension), appends a JSON receipt
IN=$1; T=$2; TAG=$3
cd /tmp/hodge/results/eight-plane-configuration/wave2-excl
OUT=${IN%.ms}.lm
START=$(date +%s)
msolve -v 1 -g 1 -t $T -f $IN -o $OUT > ${IN%.ms}.msolve.log 2>&1
RC=$?
END=$(date +%s)
DIMQ=$(grep -m1 -i "Dimension of quotient" ${IN%.ms}.msolve.log | sed 's/.*: *//')
NLM=$(grep -v '^#' $OUT | tr ',' '\n' | grep -c 'w')
echo "{\"exp\":\"b-hess-msolve\",\"input\":\"$IN\",\"rc\":$RC,\"quotient_dim\":\"$DIMQ\",\"n_leading_monomials\":$NLM,\"seconds\":$((END-START)),\"threads\":$T,\"job\":\"$TAG\",\"ts\":\"$(date +%Y-%m-%dT%H:%M:%S)\"}" >> wave2-excl-remote-$TAG.jsonl
echo "done $IN rc=$RC quotient_dim=$DIMQ n_lm=$NLM in $((END-START))s"
