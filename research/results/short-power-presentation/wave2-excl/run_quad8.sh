cd /tmp/hodge/results/short-power-presentation/wave2-excl
while kill -0 $(cat /tmp/hodge/jobs/excl-spp-quad7.pid) 2>/dev/null; do sleep 60; done
while [ ! -f e1kernel_d8_p0.json ] || [ ! -f e1kernel_d8_p1.json ]; do sleep 60; done
for m in e o; do
  python3 w2x_quad.py 8 $m > quad_d8_$m.out 2>&1
  ( time msolve -v 1 -t 6 -g 2 -f quad_d8_$m.ms -o quad_d8_$m.msolve.out ) > quad_d8_$m.msolve.log 2>&1
  python3 w2x_gbcheck.py quad_d8_$m.ms quad_d8_$m.msolve.out > quad_d8_$m.gbcheck.json 2>&1
done
