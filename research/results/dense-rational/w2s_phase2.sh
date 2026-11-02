#!/bin/bash
# phase 2 (after prep): Hessian relaxations d=7,8,9 (pure + mixed, real + complex), structured candidates, then d=8 E2 varpro
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
export W2S_THREADS=1
python3 w2s_hessrelax.py --d 7 --sector 0 --mode real --restarts 3000 --seed0 0 --procs 8 --name hr7e --max-minutes 40 > hr7e.out 2> hr7e.err
python3 w2s_hessrelax.py --d 7 --sector 1 --mode real --restarts 3000 --seed0 0 --procs 8 --name hr7o --max-minutes 40 > hr7o.out 2> hr7o.err
python3 w2s_hessrelax.py --d 7 --sector 0 --mode complex --restarts 2000 --seed0 50000 --procs 8 --name hr7ec --max-minutes 40 > hr7ec.out 2> hr7ec.err
python3 w2s_hessrelax.py --d 7 --sector 1 --mode complex --restarts 2000 --seed0 50000 --procs 8 --name hr7oc --max-minutes 40 > hr7oc.out 2> hr7oc.err
python3 w2s_hessrelax.py --d 7 --sector mixed --mode complex --restarts 2000 --seed0 70000 --procs 8 --name hr7mc --max-minutes 60 > hr7mc.out 2> hr7mc.err
python3 w2s_hessrelax.py --d 8 --sector 0 --mode real --restarts 3000 --seed0 0 --procs 8 --name hr8e --max-minutes 60 > hr8e.out 2> hr8e.err
python3 w2s_hessrelax.py --d 8 --sector 1 --mode real --restarts 3000 --seed0 0 --procs 8 --name hr8o --max-minutes 60 > hr8o.out 2> hr8o.err
python3 w2s_hessrelax.py --d 8 --sector 0 --mode complex --restarts 1500 --seed0 50000 --procs 8 --name hr8ec --max-minutes 60 > hr8ec.out 2> hr8ec.err
python3 w2s_hessrelax.py --d 8 --sector 1 --mode complex --restarts 1500 --seed0 50000 --procs 8 --name hr8oc --max-minutes 60 > hr8oc.out 2> hr8oc.err
python3 w2s_hessrelax.py --d 9 --sector 0 --mode real --restarts 1500 --seed0 0 --procs 8 --name hr9e --max-minutes 60 > hr9e.out 2> hr9e.err
python3 w2s_hessrelax.py --d 9 --sector 1 --mode real --restarts 1500 --seed0 0 --procs 8 --name hr9o --max-minutes 60 > hr9o.out 2> hr9o.err
W2S_THREADS=8 python3 w2s_structured.py > w2s_structured.out 2> w2s_structured.err
echo done > phase2.done
