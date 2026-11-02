#!/bin/bash
# local light batch B: Koszul varpro (sparse reduction, graceful timeouts) over the complete d=7 families and d=8 sectors
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
for s in 00 01 10 11; do timeout 170 python3 kvp.py --d 7 --sec $s --family J --restarts 12 --seed 17 --tag localB 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|SIGTERM" | cut -c1-250; done
for s in e o; do timeout 170 python3 kvp.py --d 7 --sec $s --family J --restarts 8 --seed 18 --tag localB 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|SIGTERM" | cut -c1-250; done
for s in 01 10 11; do timeout 170 python3 kvp.py --d 8 --sec $s --family J --restarts 8 --seed 19 --tag localB 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|SIGTERM" | cut -c1-250; done
timeout 170 python3 kvp.py --d 8 --sec 00 --family JK --restarts 8 --seed 20 --tag localB 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|SIGTERM" | cut -c1-250
echo LOCAL_BATCH_B_DONE
