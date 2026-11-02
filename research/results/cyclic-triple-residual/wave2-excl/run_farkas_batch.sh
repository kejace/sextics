#!/bin/bash
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
D8=$(python3 -c "import json; print(json.load(open('/tmp/hodge/results/cyclic-triple-residual/kernel_d8_0_4_0.json'))['D'])")
EXF=$(python3 -c "
import sys; sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
F,A,q,rec=core.load_instance('cyclic-triple-residual')
D={}
for i in range(3): D=core.add(D, core.mul({tuple(int(j==i) for j in range(6)):1}, core.pderiv(F,i)))
print(core.to_str(D))")
python3 a2_farkas.py fqlx_cone2 "X1**6 + X2**6 + Y1**6 + Y2**6"
python3 a2_farkas.py fermat_join_cone "X0**6 + X1**6 + X2**6"
python3 a2_farkas.py E_X_F "$EXF"
python3 a2_farkas.py D8 "$D8"
python3 a2_farkas.py control_F_illegal "$(python3 -c "
import sys; sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
F,A,q,rec=core.load_instance('cyclic-triple-residual'); print(core.to_str(F))")"
echo ALLDONE
