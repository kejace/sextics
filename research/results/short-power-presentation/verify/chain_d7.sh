cd /tmp/hodge/results/short-power-presentation/verify
until ! pgrep -f d7_exact.py >/dev/null && [ -f d7_exact.out ]; do sleep 20; done
if [ -f d7_N_exact.json ]; then python3 d7_e2.py > d7_e2.out 2>&1; fi
until [ "$(cut -d. -f1 /proc/loadavg)" -lt 140 ]; do sleep 60; done
python3 fullsys_v.py /tmp/hodge/results/short-power-presentation/e1vec_d7_p0.json 1000000007 D7plus > fullsys_D7.out 2>&1
