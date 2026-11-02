#!/bin/bash
# submit queue 1 (waits for remote capacity); run in background
cd /home/kejace/devel/sextics-work
L=eight-plane-configuration
sub() { name=$1; cpus=$2; shift 2; echo "[$(date +%H:%M:%S)] submitting $name"; python3 -m hodge.rjob submit --lane $L --name $name --cpus $cpus --wait 14400 -- "$@" 2>&1 | tail -2; sleep 15; }
sub search-$L-hessA 8 "python3 -u w2_hess.py --degrees 6,7 --pars 0,1 --macaulay-max 3 --macaulay-maxcols 20000 --msolve 0 --tag hessA > hessA.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 7 --par 0 --precompute 1 --stage1 1 --threads 8 > prec-d7.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 7 --par 0 --restarts 2000 --seed 1 --nproc 8 --stage1-only 1 --stage1-nfev 150 --tag stage1-d7-p0 > s1-d7-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 7 --par 1 --restarts 2000 --seed 1 --nproc 8 --stage1-only 1 --stage1-nfev 150 --tag stage1-d7-p1 > s1-d7-p1.out 2>&1"
sub search-$L-prepk 8 "python3 -u w2_k10prep.py --d 10 --par 0 --threads 8 --K 1044 > k10prep_p0.out 2>&1; python3 -u w2_k10prep.py --d 10 --par 1 --threads 8 --K 1008 > k10prep_p1.out 2>&1; OPENBLAS_NUM_THREADS=1 python3 -u w2_kernel.py --degrees 11,12 --pars 0,1 --tag kernel > kernel.out 2>&1"
sub search-$L-d6d7 8 "python3 -u w2_jvarpro.py --family jacobian --d 6 --par 0 --restarts 600 --seed 1 --nproc 8 --max-nfev 120 > s-d6-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 6 --par 1 --restarts 600 --seed 1 --nproc 8 --max-nfev 120 > s-d6-p1.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 7 --par 0 --restarts 300 --seed 1 --nproc 8 --max-nfev 150 > s-d7-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 7 --par 1 --restarts 300 --seed 1 --nproc 8 --max-nfev 150 > s-d7-p1.out 2>&1"
sub search-$L-hessB 2 "python3 -u w2_hess.py --degrees 8,9 --pars 0,1 --macaulay-max 2 --msolve 0 --tag hessB > hessB.out 2>&1"
echo "[$(date +%H:%M:%S)] queue 1 done"
