#!/bin/bash
cd /tmp/hodge/results/polarized-k3-remainder
f='X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6'
fV='Y0**6 + 5*Y0**5*Y1/7 + 5*Y0**2*Y1**3*Y2/7 + 15*Y0*Y2**5/7 + Y1**6 + 10*Y1**5*Y2/7 + Y2**6'
/usr/bin/python3 pfsolve.py U0 "X0" 2>&1 | grep -v "NOTE"
/usr/bin/python3 pfsolve.py fU "$f" 2>&1 | grep -v "NOTE"
/usr/bin/python3 pfsolve.py fU-minus-fV "($f) - ($fV)" 2>&1 | grep -v "NOTE"
/usr/bin/python3 pfsolve.py fU-times-U0 "($f)*X0" 2>&1 | grep -v "NOTE"
/usr/bin/python3 pfsolve.py fU-times-U0-minus-V0 "($f)*(X0-Y0)" 2>&1 | grep -v "NOTE"
/usr/bin/python3 pfsolve.py U0-minus-V0 "X0-Y0" 2>&1 | grep -v "NOTE"
echo PFBATCH1-DONE
