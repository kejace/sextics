"""(b) degree-7 exclusion, own replication.  Build the L-relaxation quadrics exactly over Q via Singular normal
forms modulo L = (G) + J(G)^2 + (R_ijk) (so no large exact nullspaces), gauge-fixed (drop the 6 coefficients of
z_5 z_j in a_5; the relaxation condition is gauge invariant: B_{D+tF}(T(a+tz/6)) = B_D(Ta) mod (F)), and emit
(1) quadrics.json (exact integer quadrics in the 120 gauge-fixed parameters), (2) msolve inputs mod p."""
import sys, time, json, subprocess, os
sys.path.insert(0, '.')
from relax import *
from logutil import log
d = 7
t0 = time.time()
Rx = Relaxation(d)
gauge = {pid for pid, par in enumerate(Rx.params) if par[1] == 5 and par[2][5] >= 1}
keep = [pid for pid in range(Rx.npar) if pid not in gauge]
newidx = {pid: i for i, pid in enumerate(keep)}
print('params', Rx.npar, 'gauge-fixed', len(gauge), 'remaining', len(keep))
keys = [k for k in Rx.quad_keys() if k[0] not in gauge and k[1] not in gauge]
print('quadratic monomials after gauge:', len(keys))
# T_k polys (integer) in z
Tpolys = {k: {} for k in keys}
for e, v in Rx.E2.items():
    for k, c in v.items():
        if k in Tpolys: Tpolys[k][e] = c
names = ['U0', 'U1', 'U2', 'V0', 'V1', 'V2']
def spoly(P):
    terms = []
    for e, c in P.items():
        mon = '*'.join(f'{names[i]}^{e[i]}' for i in range(NV) if e[i]) or '1'
        terms.append(f'({c})*{mon}')
    return '+'.join(terms) if terms else '0'
lines = ['ring r = 0,(U0,U1,U2,V0,V1,V2),dp;', 'option(redSB);', f'poly G = {spoly(Rx.G)};', 'ideal J = jacob(G);', 'ideal L = G;',
         'int i; int j; int k;', 'for (i=1;i<=6;i++){ for (j=i;j<=6;j++){ L = L + J[i]*J[j]; } }',
         'matrix H = jacob(J);',
         'for (i=1;i<=6;i++){ for (j=1;j<=6;j++){ for (k=j+1;k<=6;k++){ L = L + (J[k]*H[i,j] - J[j]*H[i,k]); } } }',
         'L = std(L);', 'int ne = size(L); print("GB size " + string(ne));',
         'intvec hv = hilb(L,1); print("hilbert numerator computed");',
         f'ideal KB = kbase(L,{2*d}); print("std monomials in degree {2*d}: " + string(size(KB)));',
         f'write(":w {W}/deg7_kbase.txt", string(KB));'.replace('{W}', os.getcwd()),
         f'link fh = ":w {os.getcwd()}/deg7_nf.txt";', 'poly t; poly n;']
for k in keys:
    lines.append(f't = {spoly(Tpolys[k])}; n = reduce(t, L); write(fh, "{k[0]},{k[1]}:" + string(n));')
lines.append('close(fh); print("done"); quit;')
open('deg7_nf.sing', 'w').write('\n'.join(lines) + '\n')
print('wrote deg7_nf.sing', os.path.getsize('deg7_nf.sing') // 1024, 'KB', f'{time.time()-t0:.1f}s', flush=True)
t1 = time.time()
r = subprocess.run(['Singular', '-q', 'deg7_nf.sing'], capture_output=True, text=True)
print(r.stdout[-2000:], r.stderr[-2000:], f'Singular {time.time()-t1:.1f}s', flush=True)
# parse normal forms
import re
from fractions import Fraction as Fr
kb = open('deg7_kbase.txt').read().replace('\n', '').split(',')
std = [s.strip() for s in kb if s.strip()]
print('standard monomials:', len(std))
def parse_poly(s):
    """Singular poly string over Q with vars U0..V2 -> {exp: Fraction}"""
    out = {}
    s = s.strip()
    if s == '0' or not s: return out
    for term in re.split(r'(?=[+-])', s.replace(' ', '')):
        if not term: continue
        sign = -1 if term.startswith('-') else 1
        term = term.lstrip('+-')
        coef = Fr(1); e = [0] * NV
        for fac in term.split('*'):
            m = re.fullmatch(r'([UV][012])(?:\^(\d+))?', fac)
            if m:
                e[names.index(m.group(1))] += int(m.group(2) or 1)
            else:
                coef *= Fr(fac)
        out[tuple(e)] = out.get(tuple(e), 0) + sign * coef
    return out
quad = {}   # std monomial -> {newkey: Fraction}
nnz = 0
for line in open('deg7_nf.txt'):
    line = line.strip()
    if not line: continue
    kk, s = line.split(':', 1)
    p1, p2 = map(int, kk.split(','))
    nk = (newidx[p1], newidx[p2])
    for e, c in parse_poly(s).items():
        if c: quad.setdefault(e, {})[nk] = c; nnz += 1
print('quadrics (nonzero std-monomial rows):', len(quad), 'nnz', nnz, flush=True)
# clear denominators per quadric
from math import gcd
qint = []
for e, row in quad.items():
    den = 1
    for c in row.values(): den = den * c.denominator // gcd(den, c.denominator)
    g = 0
    rowi = {k: int(c * den) for k, c in row.items()}
    for v in rowi.values(): g = gcd(g, abs(v))
    qint.append({f'{k[0]},{k[1]}': v // g for k, v in rowi.items()})
json.dump(dict(nvars=len(keep), params=[Rx.params[pid] for pid in keep], quadrics=qint), open('deg7_quadrics.json', 'w'))
# exact rank of the quadric span
rows = {}
for qd in qint:
    for k in qd: rows.setdefault(k, len(rows))
Mq = fmpz_mat(len(qint), len(rows))
for i, qd in enumerate(qint):
    for k, v in qd.items(): Mq[i, rows[k]] = v
rk = Mq.rank()
print('rank of quadric span over Q:', rk, 'of', len(qint), 'quadrics in', len(rows), 'monomials', flush=True)
log(exp='a6_deg7_build', nvars=len(keep), n_quadrics=len(qint), rank_quadrics_Q=rk, n_monomials=len(rows), std_monomials_deg14=len(std),
    seconds=round(time.time() - t0, 1), rigor='exact over Q (Singular GB of L over Q + normal forms; integer quadrics)')
