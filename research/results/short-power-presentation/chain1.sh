cd /tmp/hodge/results/short-power-presentation
python3 test_D.py e1vec_d7_p0.json --image --family "E1-kernel d=7 even (D7)" > testD7.out 2>&1
python3 e2_mod_J.py e1vec_d7_p0.json 7 0 > e2modJ_d7_p0.build.out 2>&1
python3 e1vec.py 7 1 > e1vec_7o.out 2>&1
