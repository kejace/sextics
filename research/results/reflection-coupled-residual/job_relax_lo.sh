#!/bin/bash
# low-degree Hessian-relaxation chain: exact msolve at d=7 (all sectors + mixed), LM at d=7, d=8 (sectors, mixed, affine K8)
set -x
for s in 00 01 10 11 e o; do timeout 7200 python3 rcr2.py msolve --d 7 --sec $s --run --threads 8 --seed 5; done
for s in 00 01 10 11 e o; do python3 rcr2.py relax --d 7 --sec $s --restarts 2000 --seed 11; done
python3 rcr2.py relax --d 8 --sec 00 --affine kernel_d8_X03Y03.json:0 --restarts 3000 --seed 21
python3 rcr2.py msolvef --d 8 --sec 00 --family JK --run --threads 8 --timeout 7200
python3 rcr2.py msolvef --d 8 --sec 11 --family JK --run --threads 8 --timeout 3600
python3 rcr2.py msolvef --d 8 --sec 01 --family JK --run --threads 8 --timeout 3600
python3 rcr2.py msolvef --d 8 --sec 10 --family JK --run --threads 8 --timeout 3600
for s in 00 01 10 11; do python3 rcr2.py relax --d 8 --sec $s --restarts 1500 --seed 31; done
for s in e o; do python3 rcr2.py relax --d 8 --sec $s --restarts 1000 --seed 41; done
for s in 00 01 10 11 e o; do timeout 7200 python3 rcr2.py msolve --d 8 --sec $s --run --threads 8 --seed 5; done
echo RELAX_LO_DONE
