import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
from flint import nmod_mat
import numpy as np
p = core.PRIMES[0]
F, A, q, rec = load()
Fp = core.reduce_mod(F, p)
# (1) is B = (L+M)^6 + (L-M)^6 in J(F)_6^+ ?  (equivalently e(X)+e(Y) in J)
B = core.add(core.power(core.add(L, M), 6), core.power(core.sub(L, M), 6))
Bp = core.reduce_mod(B, p)
Fi = [core.pderiv(Fp, i, p) for i in range(6)]
cols = []
for i in range(6):
    pi = 0 if i < 3 else 1
    for m in core.monomials(1, (0 - pi) % 2):
        cols.append(core.shift(Fi[i], m, 1, p))
rows = {}
for P in cols + [Bp]:
    for e in P: rows.setdefault(e, len(rows))
def rank_of(polys):
    R = len(rows); C = len(polys); flat = [0]*(R*C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e]*C+j] = v
    return nmod_mat(R, C, flat, p).rank()
rJ = rank_of(cols); rB = rank_of(cols + [Bp])
print('(1) dim J_6^+ =', rJ, '; B in J_6^+ ?', rJ == rB)
eXY = core.add(eX, eY); eXYp = core.reduce_mod(eXY, p)
print('    e(X)+e(Y) in J_6^+ ?', rank_of(cols + [eXYp]) == rJ)
# (2) Fermat P^3's in P^7 meeting Lambda: forms l_1..l_6 = coordinates, l_7 = k(L+M), l_8 = k(L-M), k = c^(1/6)
k = (5.0 / 7**7) ** (1/6)
Lv = np.array([1,2,3,0,0,0], float); Mv = np.array([0,0,0,1,3,5], float)
forms = [np.eye(6)[i] for i in range(6)] + [k*(Lv+Mv), k*(Lv-Mv)]
zetas = [np.exp(1j*np.pi*(2*j+1)/6) for j in range(6)]   # zeta^6 = -1
def pairings(items):
    if not items: yield (); return
    a = items[0]
    for b in items[1:]:
        for r in pairings([i for i in items if i not in (a,b)]): yield ((a,b),)+r
from collections import Counter
cnt = Counter(); examples = {}
t0 = time.time()
for pr in pairings(list(range(8))):
    for zs in itertools.product(range(6), repeat=4):
        rowsM = [forms[a] - zetas[z]*forms[b] for (a,b), z in zip(pr, zs)]
        r = np.linalg.matrix_rank(np.array(rowsM), tol=1e-9)
        dim = 5 - r  # projective dimension of Lambda ∩ P^3
        cnt[dim] += 1
        if dim >= 2 and dim not in examples: examples[dim] = (pr, zs)
print('(2) Fermat P^3 ∩ Lambda projective dimensions:', dict(cnt), f'({time.time()-t0:.0f}s)')
print('    examples of dim>=2:', examples)
# (3) hyperplanes V = zeta U with zeta^6=-1: F restricted should equal e(X)+e(Y)
import sympy as s
Z = s.symbols('X0 X1 X2 Y0 Y1 Y2')
Ls = Z[0]+2*Z[1]+3*Z[2]; Ms = Z[3]+3*Z[4]+5*Z[5]
zeta = s.exp(s.I*s.pi/6)
# U = L+M, V = L-M, V = zeta U  <=>  M = (1-zeta)/(1+zeta) L ;  check (1+t)^6 + (1-t)^6 = 0 for t = (1-zeta)/(1+zeta)
t = s.simplify((1-zeta)/(1+zeta))
print('(3) t = M/L on H_zeta:', s.nsimplify(s.simplify(t)), ' (1+t)^6+(1-t)^6 =', s.simplify(s.expand((1+t)**6+(1-t)**6)))
