#!/bin/bash
cd /tmp/hodge/results/klein-coupled-character/wave2-excl
for par in 0 1; do for k in 0 1 2 3 4 5 6; do
  nice -n 10 timeout 10800 python3 -u w2_radical8.py $par $k 3 > local_radical8_${par}_${k}.out 2>&1
done; done
echo ALL8_DONE >> local_radical8_1_6.out
