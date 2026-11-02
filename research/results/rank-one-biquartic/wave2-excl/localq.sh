#!/bin/bash
# sequential local queue (one heavy process at a time); each step appends to localq.log
cd /tmp/hodge/results/rank-one-biquartic/wave2-excl
run() { echo "[$(date +%H:%M:%S)] START $*" >> localq.log; "$@"; echo "[$(date +%H:%M:%S)] END rc=$? $*" >> localq.log; }
run timeout 1800 python3 hess_pipeline.py 6 3 --sparse --tag hpipe_d6 > hpipe_d6.out 2>&1
run timeout 3600 python3 k7exact.py 7 > k7exact_d7.out 2>&1
run timeout 10800 python3 hess_pipeline.py 7 3 --sparse --nonj nonj_d7.json --tag hpipe_d7 > hpipe_d7.out 2>&1
run timeout 3000 Singular -q quotJA_p65521.sing > quotJA_p65521.out 2>&1
run timeout 2400 python3 verify_e1kernel.py 1000000007 > verify_e1kernel.out 2>&1
