#!/bin/bash
# rigorous degree-7 tau-pure Hessian relaxation: msolve on the wave-1 systems (483 quadrics in 60 vars + random affine slice)
cd /tmp/hodge/results/dense-rational
export W2S_REMOTE=1
for tag in pd0_slice:polar_d7_p1073741827.ms pd1_slice:polar_d7_pd1_slice_p1073741827.ms pd0_homog:polar_d7_pd0_homog_p1073741827.ms pd1_homog:polar_d7_pd1_homog_p1073741827.ms; do
  name=${tag%%:*}; file=${tag#*:}
  echo "=== $name $file $(date)" >> w2s_ms7.log
  /usr/bin/time -v timeout 36000 msolve -v 2 -t 8 -f $file -o w2s_ms7_$name.out >> w2s_ms7_$name.msolve.log 2>&1
  echo "exit $? $(date)" >> w2s_ms7.log
  echo "--- output head:" >> w2s_ms7.log; head -c 400 w2s_ms7_$name.out >> w2s_ms7.log; echo >> w2s_ms7.log
done
echo done >> w2s_ms7.log
