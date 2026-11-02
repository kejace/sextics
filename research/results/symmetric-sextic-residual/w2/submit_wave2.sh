#!/bin/bash
# second batch (submitted when slots free up): jac w-space, S22 probes, structured families, d=12-13 invariant
cd /home/kejace/devel/sextics-work
L=symmetric-sextic-residual
python3 -m hodge.rjob submit --lane $L --name search-$L-jac1 --cpus 8 --wait 14400 -- python3 w2_varpro.py --job jac1 --sectors jac:5:0,jac:5:1,jac:6:0,jac:6:1,jac:7:0,jac:7:1 --restarts 800 --per-dim 40 --workers 8 --maxiter 400 --seed 3000
sleep 20
python3 -m hodge.rjob submit --lane $L --name search-$L-jac2 --cpus 8 --wait 14400 -- python3 w2_varpro.py --job jac2 --sectors jac:8:0,jac:8:1 --restarts 160 --workers 8 --maxiter 400 --seed 4000
sleep 20
python3 -m hodge.rjob submit --lane $L --name search-$L-struct --cpus 8 --wait 14400 -- python3 w2_struct.py --job struct --families ternX,ternY,uuv,evv,line_tt,line_sX,line_sY --chis triv,sgnX,sgnY,sgnsgn --dmin 7 --dmax 12 --restarts 150 --workers 8
sleep 20
python3 -m hodge.rjob submit --lane $L --name search-$L-vp3 --cpus 8 --wait 14400 -- python3 w2_varpro.py --job vp3 --sectors triv:12:0,triv:12:1,sgnX:8:1,sgnY:8:1,sgnX:9:0,sgnY:9:1,sgnsgn:10:1,sgnX:10:0,sgnX:10:1 --restarts 150 --workers 8 --maxiter 400 --seed 5000
echo SUBMIT-QUEUE2-DONE
