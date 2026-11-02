#!/bin/bash
# Exact fixed-D tests (Koszul normalization) of the thread's degree-8 relaxation family D = f(U) h, h quadric monomial,
# for both targets; plus served-space K3 projections of the (mu6-homogeneous) members.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
python3 - <<'PY' 2>&1 | grep -v "engine v2"
import sys, json, itertools, subprocess
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder/wave2-excl')
from fixedD_koszul import run, U, V, z, mul, scale
from fractions import Fraction as Fr
from hodge import core
names = ['U0','U1','U2','V0','V1','V2']
quads = [(i, j) for i in range(6) for j in range(i, 6)]
for (i, j) in quads:
    h = mul(z(i), z(j)); label = f'fU*{names[i]}{names[j]}'
    a = [scale(mul(U[k], h), Fr(1, 42)) for k in range(3)] + [{}, {}, {}]
    for kind in ('full', 'k3'):
        try:
            res = run(label + f'-{kind}', 8, a, kind, 65521, False, verbose=False)
            print(json.dumps({k: v for k, v in res.items() if k != 'blocks'}), flush=True)
            with open('fixedD-koszul-fUh-d8.jsonl', 'a') as fh: fh.write(json.dumps(res) + '\n')
        except Exception as ex:
            print('ERROR', label, kind, repr(ex), flush=True)
PY
echo "=== servedK3 for f(U)*h"; date
F7="X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6"
for h in "X0**2" "X0*X1" "X0*Y0" "Y0**2" "X0*Y1" "Y0*Y1" "X2**2" "Y2**2" "X1*Y2"; do
  python3 servedK3.py "fU*($h)" "($F7)*($h)" --out=servedK3-fUh-d8.jsonl 2>&1 | grep -v "engine v2"
done
echo "=== DONE"; date
