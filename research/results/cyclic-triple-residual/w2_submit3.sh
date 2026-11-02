#!/bin/bash
cd /home/kejace/devel/sextics-work
L=cyclic-triple-residual
submit() { name=$1; shift; echo "[$(date +%T)] submitting $name"; timeout 7800 python3 -m hodge.rjob submit --lane $L --name search-$L-$name --cpus 8 --wait 7200 -- "$@" 2>&1 | grep -v "engine v2"; echo "[$(date +%T)] $name exit ${PIPESTATUS[0]}"; }
submit jac bash w2_jacjob.sh
submit sec8 python3 w2_sector.py 8 1000 1000 0 600 8 ">" w2_sector8.out "2>" w2_sector8.err
submit sec9 python3 w2_sector.py 9 600 1000 0 900 8 ">" w2_sector9.out "2>" w2_sector9.err
echo "[$(date +%T)] all done"
