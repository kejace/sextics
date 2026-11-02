import sys, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core, fixed_d
F, A, q, rec = core.load_instance('cycle-q2')
Gw = eq.symmetry_group(F); gens = eq.generators(Gw)
print('group order', len(Gw), 'gens', gens.tolist())
kA = eq.key_of_poly(A, gens)
allT = eq.monomial_exps(12); kT = eq.keys_of(allT, gens)
print('S_12 sector size (my keys):', int((kT == kA).sum()))
S = eq.EqSystem(F, A, q, core.parse('1'), p=eq.PRIME, with_target=True, verbose=True)
print('EqSystem nT', S.nT, 'gens', S.gens.tolist(), 'kA', S.kA, 'order', S.order)
# check: are all monomials in the EqSystem target sector of the same (Z/6)^6-weight pattern as A under the FULL group?
Tm = np.array([c[2] for c in S.col_meta if c[0] == 'A'], dtype=np.int64)
wA = (np.array([[2,2,3,2,2,1]]) @ Gw.T) % 6
wT = (Tm @ Gw.T) % 6
print('all target monomials have A weights under full group:', bool(np.all(wT == wA)))
# validate EqSystem consistency against no-symmetry solver on a small case
D = core.parse('X0**3 + X1**3 + X2**3')
t0 = time.time()
r = eq.test_D_fast(F, A, q, D, 'validation', note='cubic cone, sector solver', verbose=True)
print('sector:', {k: v for k, v in r.items() if not k.startswith('_')}, time.time() - t0)
t0 = time.time()
out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=False, want_exact=False)
print('no-symmetry:', {k: v for k, v in out.items() if k in ('consistent', 'rank', 'nullity', 'rows', 'cols')}, time.time() - t0)
