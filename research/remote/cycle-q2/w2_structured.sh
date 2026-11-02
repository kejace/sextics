#!/bin/bash
# structured geometric families, all sectors, LM varpro; logs to wave2-search.log.jsonl
cd /tmp/hodge/results/cycle-q2
python3 w2_hess.py d6 400 > w2_hess_d6.out 2>&1
for d in 7 8 9 10; do
  for s in 0,1,2,3 0,1,2,4 0,1,2,5 0,1,3,4 0,1,3,5 0,1,4,5 0,2,3,4 0,2,3,5 0,2,4,5 0,3,4,5 1,2,3,4 1,2,3,5 1,2,4,5 1,3,4,5 2,3,4,5; do
    python3 w2_search.py --d $d --family cone4:$s --restarts 40 --batch 20 --maxiter 200 --nproc 8 --tag cone4_d${d}_${s//,/} >> w2_structured.out 2>&1
  done
  for L in 0,3 1,4 2,5 0,1 0,2 1,2 3,4 3,5 4,5 0,4 0,5 1,3 1,5 2,3 2,4; do
    python3 w2_search.py --d $d --family mult4:$L --restarts 40 --batch 20 --maxiter 200 --nproc 8 --tag mult4_d${d}_${L//,/} >> w2_structured.out 2>&1
  done
  for t in 0,1,2 3,4,5 0,1,3 0,1,4 0,1,5 0,2,3 0,2,4 0,2,5 1,2,3 1,2,4 1,2,5 0,3,4 0,3,5 0,4,5 1,3,4 1,3,5 1,4,5 2,3,4 2,3,5 2,4,5; do
    python3 w2_search.py --d $d --family cone3:$t --restarts 30 --batch 15 --maxiter 200 --nproc 8 --tag cone3_d${d}_${t//,/} >> w2_structured.out 2>&1
  done
done
echo STRUCTURED-DONE >> w2_structured.out
