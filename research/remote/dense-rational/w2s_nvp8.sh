#!/bin/bash
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
W2S_THREADS=1 python3 w2s_nvp.py --d 8 --pD 0 --kind jac --restarts 2000 --seed0 0 --maxiter 500 --procs 4 --name nvp8e --max-minutes 420 --gradcheck > nvp8e.out 2> nvp8e.err &
W2S_THREADS=1 python3 w2s_nvp.py --d 8 --pD 1 --kind jac --restarts 2000 --seed0 100000 --maxiter 500 --procs 4 --name nvp8o --max-minutes 420 --gradcheck > nvp8o.out 2> nvp8o.err &
wait
echo done > nvp8.done
