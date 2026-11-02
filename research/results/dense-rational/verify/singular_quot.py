"""Independent tool check with Singular: Hilbert functions of J, (J:A), and parity split, mod 32003.
Writes a Singular script and runs it; parses kbase monomials by parity for degrees 1..12."""
import sys, subprocess, time, re, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
p = int(sys.argv[1]) if len(sys.argv) > 1 else 32003
F, A, q = parse_rec()
Fp = reduce_mod(F, p)
def sing(P):
    return '+'.join(f'{c}*' + '*'.join(f'{VARS[i]}^{e[i]}' for i in range(6) if e[i]) for e, c in P.items())
Fi = [pderiv(Fp, i, p) for i in range(6)]
script = f"""
ring r = {p},(X0,X1,X2,Y0,Y1,Y2),dp;
poly F = {sing(Fp)};
ideal J = jacob(F);
ideal sJ = std(J);
poly A = X0^3*Y0^3;
ideal Q = quotient(sJ, ideal(A));
ideal sQ = std(Q);
int d;
for (d = 0; d <= 13; d++) {{
  ideal kJ = kbase(sJ, d);
  ideal kQ = kbase(sQ, d);
  int nJe = 0; int nJo = 0; int nQe = 0; int nQo = 0; int k;
  for (k = 1; k <= size(kJ); k++) {{ if ((deg(kJ[k], intvec(0,0,0,1,1,1)) mod 2) == 0) {{ nJe++; }} else {{ nJo++; }} }}
  for (k = 1; k <= size(kQ); k++) {{ if ((deg(kQ[k], intvec(0,0,0,1,1,1)) mod 2) == 0) {{ nQe++; }} else {{ nQo++; }} }}
  printf("DEG %s R_J_even %s R_J_odd %s R_Q_even %s R_Q_odd %s kernelmodJ_even %s kernelmodJ_odd %s", d, nJe, nJo, nQe, nQo, nJe-nQe, nJo-nQo);
}}
printf("VDIM_J %s", vdim(sJ));
quit;
"""
fn = f'/tmp/hodge/results/dense-rational/verify/quot_p{p}.sing'
open(fn, 'w').write(script)
t0 = time.time()
out = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=14000)
print(out.stdout); print(out.stderr[-2000:])
open(f'/tmp/hodge/results/dense-rational/verify/quot_p{p}.out', 'w').write(out.stdout + '\n' + out.stderr)
print('seconds', round(time.time() - t0, 1))
