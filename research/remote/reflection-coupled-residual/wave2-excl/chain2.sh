#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
P=2147483587
# wait for the running (00,11) export to finish
while pgrep -f "hess.py export 7 2147483587 00 1[1]" > /dev/null; do sleep 10; done
for pr in "00 10" "00 01" "01 10 11" "00 10 11" "00 01 11" "00 01 10" "00 01 10 11"; do
  tag=$(echo $pr | tr -d ' ')
  [ -f hess_d7_W${tag}_p${P}.ms ] || nice -n 10 timeout 1200 python3 hess.py export 7 $P $pr > export_d7_W${tag}.out 2>&1
done
nice -n 10 timeout 1500 python3 hessK.py 8 $P 00 ../kernel_d8_X03Y03.json > hessK_d8.out 2>&1
# GB checks: smallest first, bounded caps
for pr in "1011" "0111" "0110" "0011" "0010" "0001"; do
  nice -n 10 timeout 1300 python3 gbcheck.py hess_d7_W${pr}_p${P}.ms 2 1200 12 >> gbcheck_d7pairs.out 2>&1
done
nice -n 10 timeout 1900 python3 gbcheck.py hessK_d8_chi00_p${P}.ms 2 1800 14 >> gbcheck_d8K.out 2>&1
for pr in "011011" "001011" "000111" "000110" "00011011"; do
  nice -n 10 timeout 1300 python3 gbcheck.py hess_d7_W${pr}_p${P}.ms 2 1200 12 >> gbcheck_d7triples.out 2>&1
done
echo CHAIN2 DONE >> chain_d7mixed.log
