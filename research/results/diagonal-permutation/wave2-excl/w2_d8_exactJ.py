"""JOB excl-dp-d8J: exact (characteristic 0) verification that A*D8 lies in J(F), via Singular over QQ.

Computes a standard basis of J(F) = (F_0..F_5) over the rationals (modStd: modular with verification, falling back to
std), then reduce(A*D8, GB) exactly.  Also checks reduce(D8, GB) != 0 (D8 not in J).  Output d8J_<job>.json; log exp='d8_exactJ'.
"""
import os, sys, time, json, subprocess, re
import w2lib as W
from hodge import core
from hodge.core import pderiv, mul
from fractions import Fraction

job = os.environ.get('W2JOB', 'd8J')
D8 = core.parse(json.load(open(f'{W.LANE}/D8_e1kernel.json'))['D'])
F, A = W.F, W.A


def sing_poly(P):
    """Render a Fraction-coefficient poly in Singular syntax."""
    terms = []
    for e, c in P.items():
        c = Fraction(c)
        mon = '*'.join(f'{v}^{k}' for v, k in zip(core.VARS, e) if k)
        s = f'({c.numerator}/{c.denominator})' if c.denominator != 1 else f'({c.numerator})'
        terms.append(s + ('*' + mon if mon else ''))
    return '+'.join(terms) if terms else '0'


path = f'{W.PKT}/d8_exactJ.sing'
with open(path, 'w') as f:
    f.write('LIB "modstd.lib";\nring r = 0, (X0,X1,X2,Y0,Y1,Y2), dp;\n')
    f.write('poly F = ' + sing_poly(F) + ';\n')
    f.write('ideal J = jacob(F);\n')
    f.write('poly A = ' + sing_poly(A) + ';\n')
    f.write('poly D8 = ' + sing_poly(D8) + ';\n')
    f.write('int t0 = timer;\nideal G = modStd(J);\nattrib(G, "isSB", 1);\nprint("TIME " + string(timer - t0));\n')
    f.write('print("SIZEGB " + string(size(G)));\n')
    f.write('print("VDIM " + string(vdim(G)));\n')
    f.write('poly r1 = reduce(A*D8, G);\nprint("NF_AD8_ZERO " + string(r1 == 0));\n')
    f.write('poly r2 = reduce(D8, G);\nprint("NF_D8_ZERO " + string(r2 == 0));\n')
    f.write('quit;\n')
t0 = time.time()
try:
    pr = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=4 * 3600)
    out = pr.stdout + pr.stderr; status = pr.returncode
except subprocess.TimeoutExpired:
    out = ''; status = 'timeout'
secs = round(time.time() - t0, 1)
open(f'{W.PKT}/d8_exactJ.out', 'w').write(out[-5000:])
m = {k: (re.search(k + r' (\S+)', out).group(1) if re.search(k + r' (\S+)', out) else None) for k in ('TIME', 'SIZEGB', 'VDIM', 'NF_AD8_ZERO', 'NF_D8_ZERO')}
res = {'status': status, 'seconds': secs, **m}
json.dump(res, open(f'{W.PKT}/d8J_{job}.json', 'w'), indent=1)
verdict = ('A*D8 in J(F) over Q (exact), D8 not in J' if m['NF_AD8_ZERO'] == '1' and m['NF_D8_ZERO'] == '0' else f'inconclusive: {m}')
W.log('d8_exactJ', D='D8', d=8, prime=0, result=verdict, seconds=secs, note=f'Singular modStd over QQ; status {status}; vdim {m["VDIM"]}')
print('DONE', verdict, flush=True)
