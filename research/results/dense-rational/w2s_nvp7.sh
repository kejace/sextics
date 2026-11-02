#!/bin/bash
# d=7 Jacobian-type full-E2 varpro, both tau-parities, 4 processes each (8 cpus total)
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
W2S_THREADS=1 python3 w2s_nvp.py --d 7 --pD 0 --kind jac --restarts 4000 --seed0 0 --maxiter 400 --procs 4 --name nvp7e --max-minutes 300 --gradcheck > nvp7e.out 2> nvp7e.err &
W2S_THREADS=1 python3 w2s_nvp.py --d 7 --pD 1 --kind jac --restarts 4000 --seed0 100000 --maxiter 400 --procs 4 --name nvp7o --max-minutes 300 --gradcheck > nvp7o.out 2> nvp7o.err &
wait
echo done > nvp7.done
