#!/bin/bash
cd /home/kejace/devel/sextics-work
L=symmetric-sextic-residual
python3 -m hodge.rjob submit --lane $L --name search-$L-kernels --cpus 2 --force -- python3 w2_kernels.py --job kernels --part chars,sub,fermat
python3 -m hodge.rjob submit --lane $L --name search-$L-vp1 --cpus 8 --force -- python3 w2_varpro.py --job vp1 --sectors triv:8:0,triv:8:1,triv:9:0,triv:9:1,triv:6:1,triv:7:0,triv:7:1 --restarts 600 --per-dim 60 --workers 8 --maxiter 400 --seed 1000
python3 -m hodge.rjob submit --lane $L --name search-$L-vp2 --cpus 8 --force -- python3 w2_varpro.py --job vp2 --sectors triv:10:0,triv:10:1,triv:11:0,triv:11:1 --restarts 300 --workers 8 --maxiter 400 --seed 2000
python3 -m hodge.rjob submit --lane $L --name search-$L-single --cpus 8 --force -- python3 w2_single.py --job single --workers 8 --part kernels,struct --dmin 8 --dmax 12
echo FORCE-SUBMIT-DONE
