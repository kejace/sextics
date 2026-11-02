#!/bin/bash
# follow-up remote jobs (run when capacity allows and the first wave of jobs has finished):
#   sector-restricted d=7 (thousands of restarts per sector, cheap), d=10 sectors (fewer restarts), second seeds for d=8/9
cd /home/kejace/devel/sextics-work
L=cyclic-triple-residual
submit() { name=$1; shift; echo "[$(date +%T)] submitting $name"; timeout 7800 python3 -m hodge.rjob submit --lane $L --name search-$L-$name --cpus 8 --wait 7200 -- "$@" 2>&1 | grep -v "engine v2"; echo "[$(date +%T)] $name exit ${PIPESTATUS[0]}"; }
submit sec7 python3 w2_sector.py 7 3000 2000 0 300 8 ">" w2_sector7.out "2>" w2_sector7.err
submit sec10 python3 w2_sector.py 10 300 600 0 1200 8 ">" w2_sector10.out "2>" w2_sector10.err
submit sec8b python3 w2_sector.py 8 1000 0 5000 600 8 ">" w2_sector8b.out "2>" w2_sector8b.err
submit sec9b python3 w2_sector.py 9 600 0 5000 900 8 ">" w2_sector9b.out "2>" w2_sector9b.err
echo "[$(date +%T)] all done"
