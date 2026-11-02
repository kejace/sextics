"""Shapes test (independent re-run of repr.py) at a new prime, reporting ranks so rigor can be judged:
system T g = -K8|rows ; if rank_p T = 168 (= #columns) then rank_Q T = 168 too and 'unsolvable mod p' is rigorous over Q."""
import sys, json, time, itertools
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
Fs, As, qv = load_rec(); F = modp(parse(Fs), p)
K8 = modp(parse(json.load(open('/tmp/hodge/results/reflection-coupled-residual/kernel_d8_X03Y03.json'))[0]), p)
Fi = [deriv(F, i, p) for i in range(6)]; cFi = [char_of(f) for f in Fi]
mons = monos_chi(8, (0, 0)) + monos_chi(8, (1, 1))   # all Y-even degree-8 monomials
idx = {m: j for j, m in enumerate(mons)}; C = len(mons)
T = [shift(Fi[i], m) for i in range(6) for m in monos(3) if chi_mul(chi(m), cFi[i])[0] == chi_mul(chi(m), cFi[i])[1]]  # g_i F_i Y-even
ng = len(T); print('|S_8^even| =', C, ' #generators of J_8^even =', ng)
Tv = []
for t in T:
    v = [0]*C
    for e, c in t.items(): v[idx[e]] = c
    Tv.append(v)
Kv = [0]*C
for e, c in K8.items(): Kv[idx[e]] = c
V = ['X0','X1','X2','Y0','Y1','Y2']
results = []
def test(name, forbid):
    rows = [j for j, m in enumerate(mons) if forbid(m)]
    if not rows: return
    M = nmod_mat(len(rows), ng + 1, p); M0 = nmod_mat(len(rows), ng, p)
    for r, j in enumerate(rows):
        for k in range(ng):
            if Tv[k][j]: M[r, k] = Tv[k][j]; M0[r, k] = Tv[k][j]
        M[r, ng] = (-Kv[j]) % p
    r0 = M0.rank(); r1 = M.rank()
    ok = (r1 == r0)
    rig = 'rigorous' if (not ok and r0 == ng) else ('mod-p only' if not ok else '-')
    print(f'{name:48s} constraints={len(rows):4d} rankT={r0:3d}/{ng} solvable={ok} transfer={rig}', flush=True)
    results.append({'shape': name, 'constraints': len(rows), 'rankT': r0, 'solvable_modp': ok, 'transfer_to_Q': rig})
test('divisible by Y2^2 (D = u*G)', lambda m: m[5] < 2)
test('divisible by Y2^4', lambda m: m[5] < 4)
test('no Y2 (D in C[x])', lambda m: m[5] > 0)
for i in range(6):
    test(f'divisible by {V[i]}', lambda m, i=i: m[i] == 0)
    test(f'divisible by {V[i]}^2', lambda m, i=i: m[i] < 2)
for i, j in itertools.combinations(range(6), 2):
    test(f'independent of {V[i]},{V[j]} (cone, vertex line)', lambda m, i=i, j=j: m[i] > 0 or m[j] > 0)
for i, j in [(0,1),(3,4),(0,3),(0,5),(3,5),(1,4)]:
    test(f'in ideal ({V[i]},{V[j]})^2', lambda m, i=i, j=j: m[i] + m[j] < 2)
    test(f'divisible by {V[i]}*{V[j]}', lambda m, i=i, j=j: m[i] == 0 or m[j] == 0)
test('Y-degree <= 2', lambda m: m[3] + m[4] + m[5] > 2)
test('Y-degree <= 4', lambda m: m[3] + m[4] + m[5] > 4)
test('X-degree <= 2', lambda m: m[0] + m[1] + m[2] > 2)
test('X-degree <= 4', lambda m: m[0] + m[1] + m[2] > 4)
test('only even exponents (poly in squares)', lambda m: any(x % 2 for x in m))
test('in ideal (X0,X1,X2)^4', lambda m: m[0]+m[1]+m[2] < 4)
test('in ideal (Y0,Y1,Y2)^4', lambda m: m[3]+m[4]+m[5] < 4)
test('in ideal (X0,Y0)^4', lambda m: m[0]+m[3] < 4)
test('in ideal (X0,Y0)^6', lambda m: m[0]+m[3] < 6)
test('in (X0,X1,Y0,Y1)^6', lambda m: m[0]+m[1]+m[3]+m[4] < 6)
test('in (X2,Y2)^2', lambda m: m[2]+m[5] < 2)
test('in (X2,Y2)^4', lambda m: m[2]+m[5] < 4)
# extra shapes the note did not test
test('independent of Y2 and Y0', lambda m: m[5] > 0 or m[3] > 0)
test('a polynomial in (X0,X1,X2,Y0*Y1,Y0^2,Y1^2,Y2^2) i.e. rho-invariant monomials only', lambda m: (m[3]+m[4]) % 2 or m[5] % 2)
json.dump(results, open(f'shapes_p{p}.json', 'w'), indent=1)
logjson(VLOG, {'check': 'shapes', 'p': p, 'n_shapes': len(results), 'any_solvable': any(r['solvable_modp'] for r in results),
               'all_rigorous': all(r['transfer_to_Q'] == 'rigorous' for r in results if not r['solvable_modp'])})
