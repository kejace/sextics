cd /tmp/hodge/results/short-power-presentation/wave2-excl
for m in e o m; do python3 w2x_quad.py 6 $m --run-singular > quad_d6_$m.out 2>&1; done
