cd /tmp/hodge/results/short-power-presentation/wave2-excl
# wait for the d=6 chain and for the d=7 kernel bases from the e1kernel ladder
while kill -0 $(cat /tmp/hodge/jobs/excl-spp-quad6.pid) 2>/dev/null; do sleep 30; done
while [ ! -f e1kernel_d7_p0.json ] || [ ! -f e1kernel_d7_p1.json ]; do sleep 30; done
for m in e o m; do python3 w2x_quad.py 7 $m --run-singular > quad_d7_$m.out 2>&1; done
