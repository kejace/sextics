#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
P=2147483587
[ -f hess_d7_W000110_p${P}.ms ] || nice -n 10 timeout 1500 python3 hess.py export 7 $P 00 01 10 > export_d7_W000110.out 2>&1
[ -f hess_d7_W00011011_p${P}.ms ] || nice -n 10 timeout 1800 python3 hess.py export 7 $P 00 01 10 11 > export_d7_W00011011.out 2>&1
[ -f hessK_d8_chi00_p${P}.ms ] || nice -n 10 timeout 1800 python3 hessK.py 8 $P 00 ../kernel_d8_X03Y03.json > hessK_d8.out 2>&1
nice -n 10 timeout 1900 python3 gbcheck.py hessK_d8_chi00_p${P}.ms 2 1800 14 >> gbcheck_d8K.out 2>&1
for pr in "0110" "0011"; do
  nice -n 10 timeout 1300 python3 gbcheck.py hess_d7_W${pr}_p${P}.ms 2 1200 12 >> gbcheck_d7pairs.out 2>&1
done
echo CHAIN4 DONE >> chain_d7mixed.log
