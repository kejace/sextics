#!/bin/bash
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
W2S_THREADS=2 python3 w2s_nvp.py --d 9 --pD 0 --kind jac --restarts 600 --seed0 0 --maxiter 600 --procs 2 --name nvp9e --max-minutes 420 --gradcheck > nvp9e.out 2> nvp9e.err &
W2S_THREADS=2 python3 w2s_nvp.py --d 9 --pD 1 --kind jac --restarts 600 --seed0 100000 --maxiter 600 --procs 2 --name nvp9o --max-minutes 420 --gradcheck > nvp9o.out 2> nvp9o.err &
wait
echo done > nvp9.done
