#!/bin/bash
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
W2S_THREADS=4 python3 w2s_nvp.py --d 10 --pD 0 --kind ker --restarts 200 --seed0 0 --maxiter 600 --procs 2 --name nvp10ker --max-minutes 420 --gradcheck > nvp10ker.out 2> nvp10ker.err
W2S_THREADS=4 python3 w2s_nvp.py --d 10 --pD 0 --kind kfull --restarts 200 --seed0 0 --maxiter 800 --procs 2 --name nvp10full --max-minutes 300 > nvp10full.out 2> nvp10full.err
echo done > nvp10.done
