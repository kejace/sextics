"""Structural candidates at degrees 9,10 (subsampled monomial multipliers) + all products of partials."""
import oplib, time, sys, itertools, random, json
from hodge import core
from hodge.core import mul, scale, pderiv, add, sub, monomials, shift
F, A, q = oplib.load()
p = oplib.P0
Fi = [pderiv(F, i) for i in range(6)]
D0 = {}
for e, c in Fi[0].items():
    e2 = list(e); e2[0] -= 1; D0[tuple(e2)] = c / 6
V = core.VARS
def mstr(e): return '*'.join(f'{V[i]}**{k}' if k > 1 else V[i] for i, k in enumerate(e) if k) or '1'
rng = random.Random(2026)
cands = []
# degree 9
for i in range(6):
    cands.append((f'D0*F_{V[i]}', 9, 'quartic-colon*partial', mul(D0, Fi[i])))
ms = monomials(4); rng.shuffle(ms)
for m in ms[:20]:
    for i in range(6):
        cands.append((f'F_{V[i]}*{mstr(m)}', 9, 'partial*monomial(sample)', shift(Fi[i], m)))
ms5 = monomials(5); rng.shuffle(ms5)
for m in ms5[:25]:
    cands.append((f'D0*{mstr(m)}', 9, 'quartic-colon*monomial(sample)', shift(D0, m)))
# degree 10
for (i, j) in itertools.combinations_with_replacement(range(6), 2):
    cands.append((f'F_{V[i]}*F_{V[j]}', 10, 'partial*partial', mul(Fi[i], Fi[j])))
ms6 = monomials(6); rng.shuffle(ms6)
for m in ms6[:20]:
    cands.append((f'D0*{mstr(m)}', 10, 'quartic-colon*monomial(sample)', shift(D0, m)))
ms5b = monomials(5); rng.shuffle(ms5b)
for m in ms5b[:15]:
    for i in range(6):
        cands.append((f'F_{V[i]}*{mstr(m)}', 10, 'partial*monomial(sample)', shift(Fi[i], m)))
# D0^2 * quadratic monomials (degree 10), D0 * Hessian minors? (degree 12, skip)
for m in monomials(2):
    cands.append((f'D0**2*{mstr(m)}', 10, 'quartic-colon-square*monomial', shift(mul(D0, D0), m)))
print('candidates:', len(cands), file=sys.stderr)
seen = set(); nsolved = 0
for name, d, fam, D in cands:
    key = tuple(sorted(core.reduce_mod(D, p).items()))
    if key in seen: continue
    seen.add(key)
    if oplib.char_of(D) is None:
        print('skip non-eigen', name, file=sys.stderr); continue
    t0 = time.time()
    ok, info = oplib.e1_reduced(F, A, q, D, p)
    if not ok:
        oplib.log(name, d, fam, False, None, info['rank'], 0, time.time()-t0, f'{name}: E1 fails'); continue
    r = oplib.reduced_solve(F, A, q, D, p)
    dt = time.time() - t0
    oplib.log(core.to_str(D) if len(D) < 40 else name, d, fam, True, r['consistent'], r['rank'], r['nullity'], dt,
              f"{name}: char={oplib.char_of(D)} rows={r['rows']} cols={r['cols']}")
    if r['consistent']:
        nsolved += 1; print('!!! CONSISTENT:', name, file=sys.stderr)
        st, data = oplib.lift_solution(F, A, q, D)
        if st == 'solved':
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open(f'witness-candidate-{name.replace("*","_")}.json', 'w'))
            print('!!! WITNESS SAVED', file=sys.stderr)
    print(f'{name:40s} d={d} char={oplib.char_of(D)} consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]} {dt:.1f}s', file=sys.stderr)
print('done; consistent count', nsolved, file=sys.stderr)
