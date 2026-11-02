import oplib, time, json, sys
from hodge import core, fixed_d
from hodge.core import mul, scale, pderiv, add, sub
F, A, q = oplib.load()
# (0) D = 1
t0 = time.time()
ok, info = fixed_d.e1_solvable(F, A, q, core.parse('1'))
st, data = fixed_d.solve_fixed_D(F, A, q, core.parse('1'))
dt = time.time() - t0
print('D=1: e1', ok, info, 'full', st, {k: v for k, v in data.items() if k != 'N' and k != 'Q'})
oplib.log('1', 0, 'sanity', ok, (st == 'solved'), data.get('rank', info['rank']), data.get('nullity', 0), dt, 'A not in J(F); class nonzero')
# quartic colon direction D0 = F_X0/(6 X0)
FX0 = pderiv(F, 0)
D0 = {}
for e, c in FX0.items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = c / 6
print('D0 =', core.to_str(D0), 'char', oplib.char_of(D0))
t0 = time.time(); ok, info = fixed_d.e1_solvable(F, A, q, D0); st, data = fixed_d.solve_fixed_D(F, A, q, D0); dt = time.time() - t0
print('D0 full:', ok, info, st, {k: v for k, v in data.items() if k not in ('N', 'Q')}, f'{dt:.1f}s')
t0 = time.time(); r = oplib.reduced_solve(F, A, q, D0, verbose=True); dt2 = time.time() - t0
print('D0 reduced:', {k: v for k, v in r.items() if k not in ('x', 'pivots', 'S')}, f'{dt2:.1f}s')
oplib.log(core.to_str(D0), 4, 'colon-quartic', ok, r['consistent'], r['rank'], r['nullity'], dt + dt2, f"full ranks {data} ; reduced rows={r['rows']} cols={r['cols']}")
# positive control: served targets of D0 in full system; pick one of char (3,0) and test reduced system
t0 = time.time()
dim, basis, S = fixed_d.target_image(F, q, D0)
print('target_image(D0): dim', dim, f'{time.time()-t0:.1f}s')
cands = [b for b in basis if oplib.char_of(b) == (3, 0)]
print('served targets with char (3,0):', len(cands))
# project a random combination onto char (3,0)?  Basis vectors need not be eigen; build the char-(3,0) projection of the image space
import flint
from flint import nmod_mat
p = oplib.P0
mons = [m for m in core.monomials(6, 1) if oplib.char(m) == (3, 0)]
idx = {m: i for i, m in enumerate(mons)}
# the image space is G-stable (D0 eigenvector), so its (3,0)-isotypic part = projection of basis vectors
rows = []
for b in basis:
    v = [0]*len(mons)
    for e, c in b.items():
        if e in idx: v[idx[e]] = c % p
    rows.append(v)
M = nmod_mat(len(rows), len(mons), [x for r_ in rows for x in r_], p)
Mr, rk = M.rref()
print('isotypic (3,0) part of served targets: dim', rk, 'of', len(mons), 'monomials')
if rk:
    Aprime = {}
    for j in range(len(mons)):
        v = int(Mr[0, j])
        if v: Aprime[mons[j]] = v
    # lift to small rationals if possible
    Al = oplib.lift_poly_small(Aprime, p)
    print('A\' (mod p, first basis vector) terms:', len(Aprime), 'lift ok' if Al else 'no small lift')
    if Al:
        print("A' =", core.to_str(Al))
        st2, d2 = oplib.lift_solution(F, Al, q, D0, verbose=False)
        print('reduced system for served target:', st2, {k: v for k, v in d2.items() if k not in ('N','Q')})
        st3, d3 = fixed_d.solve_fixed_D(F, Al, q, D0)
        print('full system for served target:', st3, {k: v for k, v in d3.items() if k not in ('N','Q')})
