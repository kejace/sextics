"""Structural candidate denominators (eigenvectors) at degrees 7..10: E1 filter + reduced fixed-D solve."""
import oplib, time, sys, itertools
from hodge import core
from hodge.core import mul, scale, pderiv, add, sub, monomials, power, shift
F, A, q = oplib.load()
p = oplib.P0
Fi = [pderiv(F, i) for i in range(6)]
Fij = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
FX0 = Fi[0]; D0 = {}
for e, c in FX0.items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = c / 6
V = core.VARS
def mono(e): return {tuple(e): 1}
def mstr(e): return '*'.join(f'{V[i]}**{k}' if k > 1 else V[i] for i, k in enumerate(e) if k) or '1'
degs = [int(a) for a in sys.argv[1:]] or [7, 8]
cands = []
for d in degs:
    # F_i * monomial
    for i in range(6):
        for m in monomials(d - 5):
            cands.append((f'F_{V[i]}*{mstr(m)}', d, 'partial*monomial', shift(Fi[i], m)))
    # D0 * monomial
    for m in monomials(d - 4):
        cands.append((f'D0*{mstr(m)}', d, 'quartic-colon*monomial', shift(D0, m)))
    if d == 8:
        # Hessian 2x2 minors F_ik F_jl - F_il F_jk
        for (i, j) in itertools.combinations(range(6), 2):
            for (k, l) in itertools.combinations_with_replacement(range(6), 2):
                Dm = sub(mul(Fij[i][k], Fij[j][l]), mul(Fij[i][l], Fij[j][k]))
                if Dm: cands.append((f'Hess[{i}{j}|{k}{l}]', d, 'hessian-2x2-minor', Dm))
        cands.append(('D0**2', d, 'quartic-colon-square', mul(D0, D0)))
    if d == 9:
        for i in range(6):
            cands.append((f'D0*F_{V[i]}', d, 'quartic-colon*partial', mul(D0, Fi[i])))
    if d == 10:
        for (i, j) in itertools.combinations_with_replacement(range(6), 2):
            cands.append((f'F_{V[i]}*F_{V[j]}', d, 'partial*partial', mul(Fi[i], Fi[j])))
    # subset-Euler derivatives E_I(F) times monomials of degree d-6
    if d >= 7:
        for r in (1, 2, 3):
            for I in itertools.combinations(range(6), r):
                EI = {}
                for i in I: EI = add(EI, shift(Fi[i], tuple(1 if j == i else 0 for j in range(6))))
                for m in monomials(d - 6):
                    cands.append((f'E_{"".join(V[i] for i in I)}(F)*{mstr(m)}', d, 'subset-euler*monomial', shift(EI, m)))
print('candidates:', len(cands), file=sys.stderr)
seen = set()
nsolved = 0
for name, d, fam, D in cands:
    key = tuple(sorted(core.reduce_mod(D, p).items()))
    if key in seen: continue
    seen.add(key)
    if oplib.char_of(D) is None:
        print('skip non-eigen', name); continue
    t0 = time.time()
    ok, info = oplib.e1_reduced(F, A, q, D, p)
    if not ok:
        oplib.log(core.to_str(D) if len(D) < 40 else name, d, fam, False, None, info['rank'], 0, time.time()-t0, f'{name}: E1 fails (reduced ranks {info["rank"]}/{info["rank_aug"]})')
        continue
    r = oplib.reduced_solve(F, A, q, D, p)
    dt = time.time() - t0
    oplib.log(core.to_str(D) if len(D) < 40 else name, d, fam, True, r['consistent'], r['rank'], r['nullity'], dt,
              f"{name}: char={oplib.char_of(D)} rows={r['rows']} cols={r['cols']}")
    if r['consistent']:
        nsolved += 1
        print('!!! CONSISTENT:', name, core.to_str(D), file=sys.stderr)
        st, data = oplib.lift_solution(F, A, q, D)
        if st == 'solved':
            import json
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open(f'witness-candidate-{name.replace("*","_").replace("/","_")}.json', 'w'))
            print('!!! WITNESS SAVED', file=sys.stderr)
    print(f'{name:40s} d={d} char={oplib.char_of(D)} consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} {dt:.1f}s', file=sys.stderr)
print('done; consistent count', nsolved, file=sys.stderr)
