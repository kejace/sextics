cd /tmp/hodge/results/short-power-presentation/verify
python3 relax6.py 0 e > relax6_e_p0.out 2>&1
python3 relax6.py 0 o > relax6_o_p0.out 2>&1
python3 relax6.py 1000000007 e > relax6_e_p1e9.out 2>&1
python3 relax6.py 1000000007 o > relax6_o_p1e9.out 2>&1
until [ "$(cut -d. -f1 /proc/loadavg)" -lt 140 ]; do sleep 60; done
python3 relax6.py 1000000007 m > relax6_m_p1e9.out 2>&1
