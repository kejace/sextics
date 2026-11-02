# remote: full Groebner bases of the d=7 relaxations (random affine charts), 16 threads, generous memory
cd /tmp/hodge/results/short-power-presentation
for par in 1 0; do
  f=relax_d7_${par}_w2_chart.ms
  echo "== msolve chart d=7 par=$par $(date)"
  /usr/bin/time -v msolve -t 16 -v 1 -f $f -o relax_d7_${par}_w2_chart.remote.msout > relax_d7_${par}_w2_chart.remote.log 2>&1
  echo "rc=$? $(date)"; head -c 300 relax_d7_${par}_w2_chart.remote.msout; echo
done
echo "== done $(date)"
