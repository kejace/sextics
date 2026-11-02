#!/bin/bash
cd /tmp/hodge/results/polarized-k3-remainder
while ! grep -q SEQDRIVE-DONE drive_low.log 2>/dev/null; do sleep 30; done
f='X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6'
fV='Y0**6 + 5*Y0**5*Y1/7 + 5*Y0**2*Y1**3*Y2/7 + 15*Y0*Y2**5/7 + Y1**6 + 10*Y1**5*Y2/7 + Y2**6'
for args in "fU|$f" "fU-minus-fV|($f) - ($fV)" "U0-minus-V0|X0-Y0" "U0U1|X0*X1" "U0V0|X0*Y0" "sumU3-sumV3|X0**3+X1**3+X2**3-Y0**3-Y1**3-Y2**3"; do
  label=${args%%|*}; expr=${args#*|}
  echo "### $(date +%T) pfsolve $label" >> pf_batch2.log
  timeout 2400 /usr/bin/python3 pfsolve.py "$label" "$expr" 2>&1 | grep -v NOTE >> pf_batch2.log; echo "EXIT=$?" >> pf_batch2.log
done
echo PFBATCH2-DONE >> pf_batch2.log
# served-space + K3 projection for key d=8 candidates
/usr/bin/python3 - << 'PY' > .img_cmds
import json, shlex
C = {c[1]: c for c in json.load(open('candidates_d8_priority.json'))}
for l in ['sumU8+sumV8', 'Y0^8+Y1^8+Y2^8', 'fU-U0sq-sym', 'M01(U)+M01(V)']:
    fam, label, expr, frame = C[l]
    print(f"timeout 3000 /usr/bin/python3 leanimage.py {shlex.quote(label)} {shlex.quote(expr)} {frame}")
PY
while read -r cmd; do echo "### $(date +%T) $cmd" | cut -c1-160 >> images.log; bash -c "$cmd" 2>&1 | grep -v NOTE >> images.log; echo "EXIT=$?" >> images.log; done < .img_cmds
echo IMAGES-DONE >> images.log
