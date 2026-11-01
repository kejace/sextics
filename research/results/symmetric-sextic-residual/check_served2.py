import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
p = P0
QX = e_k(X, 2)
rk, basis, S = served_invariant_targets(F, q, QX)
print('inv served dim', rk)
for b in basis:
    # b: dict rep-monomial -> coeff (signed); expand to invariant poly
    P = {}
    for m, v in b.items():
        P = core.add(P, core.scale(orbit_sum(m), Fraction(v)))
    print('  basis elt:', core.to_str(P)[:400])
    # check J membership of the expanded poly via e1 with D=1? use engine full target_image instead
# engine full target image for D=QX
t0=time.time()
dim, fb, S2 = fixed_d.target_image(F, q, QX)
print('engine full served dim for QX:', dim, f'({time.time()-t0:.1f}s)')
# how many of those are G-invariant? project: compute rank of invariant-sym of basis
inv = []
for b in fb:
    Pb = {e: (v if v < p//2 else v - p) for e, v in b.items()}
    inv.append(core.reduce_mod(sym(Pb), p))
Treps = orbit_reps(6, 1, canon_codes)
import numpy as np
from flint import nmod_mat
M = nmod_mat(len(inv), len(Treps), [inv[i].get(r, 0) for i in range(len(inv)) for r in Treps], p)
print('rank of G-symmetrised engine-served space:', M.rank())
# exact solve: first invariant basis element as target with D = QX
P = {}
for m, v in basis[0].items():
    P = core.add(P, core.scale(orbit_sum(m), Fraction(v)))
st, data = fixed_d.solve_fixed_D(F, P, q, QX, verbose=False)
print('engine solve_fixed_D(A=basis0, D=QX):', st, {k: v for k, v in data.items() if k not in ('N','Q')})
if st == 'solved':
    ok, info = core.verify_witness(F, P, q, QX, data['N'], data['Q'])
    print('verify:', ok, info)
    print('N =', [core.to_str(n)[:150] for n in data['N']])
    print('Q =', core.to_str(data['Q'])[:200])
# also: is basis0 served by D = 1?
st1, d1 = fixed_d.solve_fixed_D(F, P, q, core.parse('1'))
print('D=1:', st1)
