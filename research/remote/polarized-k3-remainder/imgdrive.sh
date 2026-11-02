#!/bin/bash
cd /tmp/hodge/results/polarized-k3-remainder
/usr/bin/python3 - << 'PY' > .img_cmds2
import json, shlex
C = {c[1]: c for c in json.load(open('candidates_d8_priority.json'))}
for l in ['sum_i f_i(U)V_i^3 + swap', 'U0^8+U1^8+V0^8+V1^8', 'X0^8+X1^8+X2^8', 'M01(f)(Y)', 'sumU8-sumV8', 'sum_i f_i(U)U_i^3 + swap', 'M01+M02+M12(U)+(V)']:
    fam, label, expr, frame = C[l]
    print(f"timeout 5400 /usr/bin/python3 leanimage.py {shlex.quote(label)} {shlex.quote(expr)} {frame}")
PY
while read -r cmd; do echo "### $(date +%T) $cmd" | cut -c1-160 >> images2.log; bash -c "$cmd" 2>&1 | grep -v NOTE >> images2.log; echo "EXIT=$?" >> images2.log; done < .img_cmds2
echo IMAGES2-DONE >> images2.log
