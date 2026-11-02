#!/bin/bash
# degree-8 pure-character Jacobian families (chi = 01, 10, 11; the (0,0) family with K8 is handled by hessK locally/remotely)
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
export W2LIB_SUFFIX=.r3
P=2147483587
for c in 11 10 01; do
  python3 hess.py export 8 $P $c > export_d8_W$c.out 2>&1
  python3 gbcheck.py hess_d8_W${c}_p${P}.ms 4 5400 14 >> gbcheck_d8.out 2>&1
done
python3 hessK.py 8 $P 00 ../kernel_d8_X03Y03.json > hessK_d8_remote.out 2>&1
python3 gbcheck.py hessK_d8_chi00_p${P}.ms 8 7200 14 >> gbcheck_d8.out 2>&1
echo HESS8 DONE >> job_hess8.log
