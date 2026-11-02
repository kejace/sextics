# when the d=8 relaxation files are written, stop the (hopeless) full-family msolve/Singular runs and start subfamily restrictions
cd /tmp/hodge/results/short-power-presentation
for par in 0 1; do
  while [ ! -f relax_d8_${par}_w2_chart.ms ]; do sleep 20; done
  sleep 30
  pkill -f "[w]2_relax.py 8 ${par} " ; pkill -f "[m]solve.*relax_d8_${par}_w2" ; pkill -f "[S]ingular -q relax_d8_${par}_w2"
  echo "d=8 par=$par: full-family solvers stopped $(date); starting subfamily restrictions"
  python3 w2_subfam2.py relax_d8_${par}_w2 --maxr 5 --timeout 300 --tag w2 > w2_subfam_d8_${par}.out 2>&1
  echo "d=8 par=$par subfamilies done $(date)"
done
