#!/bin/bash
# submit queue 2: waits for (a) my own running-job count < 4 and (b) remote capacity (rjob --wait)
cd /home/kejace/devel/sextics-work
L=eight-plane-configuration
myjobs() { python3 -m hodge.rjob status 2>/dev/null | grep -c "RUNNING search-$L-"; }
waitslot() { while [ "$(myjobs)" -ge 4 ]; do sleep 120; done; }
sub() { name=$1; cpus=$2; shift 2; waitslot; echo "[$(date +%H:%M:%S)] submitting $name (my running: $(myjobs))"; python3 -m hodge.rjob submit --lane $L --name $name --cpus $cpus --wait 14400 -- "$@" 2>&1 | tail -2; sleep 20; }
# d=8 / d=9: precompute the R-projector once per degree (8 threads), then searches
sub search-$L-d8 8 "python3 -u w2_jvarpro.py --family jacobian --d 8 --par 0 --precompute 1 --stage1 1 --threads 8 > prec-d8.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 8 --par 0 --restarts 160 --seed 1 --nproc 8 --max-nfev 150 --stage1 1 > s-d8-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 8 --par 1 --restarts 160 --seed 1 --nproc 8 --max-nfev 150 --stage1 1 > s-d8-p1.out 2>&1"
sub search-$L-fam10 2 "python3 -u w2_fam10.py --d 10 --pars 0,1 --tag fam10 > fam10.out 2>&1"
sub search-$L-d9 8 "python3 -u w2_jvarpro.py --family jacobian --d 9 --par 0 --precompute 1 --stage1 1 --threads 8 > prec-d9.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 9 --par 0 --restarts 400 --seed 1 --nproc 8 --stage1-only 1 --stage1-nfev 150 --tag stage1-d9-p0 > s1-d9-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 9 --par 1 --restarts 400 --seed 1 --nproc 8 --stage1-only 1 --stage1-nfev 150 --tag stage1-d9-p1 > s1-d9-p1.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 9 --par 0 --restarts 48 --seed 1 --nproc 8 --max-nfev 100 --stage1 1 > s-d9-p0.out 2>&1; python3 -u w2_jvarpro.py --family jacobian --d 9 --par 1 --restarts 48 --seed 1 --nproc 8 --max-nfev 100 --stage1 1 > s-d9-p1.out 2>&1"
sub search-$L-k10 8 "while [ ! -f w2_kernelfam_d10_p0.npz ]; do sleep 60; done; python3 -u w2_jvarpro.py --family kernel --d 10 --par 0 --npz w2_kernelfam_d10_p0.npz --restarts 24 --seed 1 --nproc 8 --max-nfev 80 > s-k10-p0.out 2>&1; while [ ! -f w2_kernelfam_d10_p1.npz ]; do sleep 60; done; python3 -u w2_jvarpro.py --family kernel --d 10 --par 1 --npz w2_kernelfam_d10_p1.npz --restarts 24 --seed 1 --nproc 8 --max-nfev 80 > s-k10-p1.out 2>&1"
echo "[$(date +%H:%M:%S)] queue 2 done"
