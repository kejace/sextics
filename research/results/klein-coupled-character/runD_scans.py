"""Exhaustive 1-parameter scans over GF(p) for small p (rational parameter values are caught up to finitely many bad primes),
plus a handful of structured candidates."""
import sys, time, random, json
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
F, A, q = load()
Fi = [core.pderiv(F, i) for i in range(6)]
X = [core.parse(v) for v in core.VARS]
def lin(polys, coeffs):
    acc = {}
    for P, c in zip(polys, coeffs):
        if c: acc = core.add(acc, core.scale(P, Fraction(c)))
    return acc
radX = lin([core.mul(X[i], Fi[i]) for i in range(3)], [1, 1, 1])   # X.grad_X F = 6K + 4M
radY = lin([core.mul(X[i], Fi[i]) for i in range(3, 6)], [1, 1, 1]) # Y.grad_Y F = 6G + 2M
print('radX terms', len(radX), 'char', sigma_weight(radX), 'radY terms', len(radY), 'char', sigma_weight(radY), flush=True)

def scan(family_name, make_D, params, p, note=''):
    """make_D(t) -> D (rational or mod p dict); scan all t in params over GF(p)."""
    t0 = time.time(); hits = []; ranks = {}
    ctx = Ctx(F, A, p)
    for t in params:
        D = make_D(t)
        if not D: continue
        ok, N0, piv = e1_lift(ctx, D, True)
        if not ok:
            ranks.setdefault('e1fail', 0); ranks['e1fail'] += 1; continue
        res = koszul_test(ctx, D, N0, True)
        key = (res['consistent'], res['rank'])
        ranks[str(key)] = ranks.get(str(key), 0) + 1
        if res['consistent']: hits.append(t)
    el = round(time.time() - t0, 1)
    log({'D': f'scan {family_name} over GF({p})', 'd': core.degree(make_D(params[1])), 'family': 'scan-1param', 'engine': 'koszul',
         'e1': None, 'consistent': bool(hits), 'rank': None, 'nullity': None, 'seconds': el,
         'note': f'{len(params)} parameter values; hits={hits[:20]}; rank histogram {ranks}; {note}'})
    print(f'SCAN {family_name} p={p}: hits={hits[:20]} ranks={ranks} [{el}s]', flush=True)
    return hits

# (i) radial family  X.grad_X F + t Y.grad_Y F  (t = 1 gives 6F: skipped by gcd; still tested, expect degenerate)
for p in [101, 103]:
    scan('radX + t*radY', lambda t: core.add(core.reduce_mod(radX, p), core.scale(core.reduce_mod(radY, p), t, p), p), list(range(p)), p,
         note='t=1 is 6F (degenerate)')
    scan('t*radX + radY', lambda t: core.add(core.scale(core.reduce_mod(radX, p), t, p), core.reduce_mod(radY, p), p), [0], p)
# (ii) Y-odd polars: F_3 + t F_4 (char 0) ; F_3 + t F_5 ; F_4 + t F_5 ; and the 2-param family F_3 + t F_4 + u F_5 over GF(101)
for p in [101]:
    scan('F_3 + t F_4', lambda t: core.add(core.reduce_mod(Fi[3], p), core.scale(core.reduce_mod(Fi[4], p), t, p), p), list(range(p)), p)
    scan('F_3 + t F_5', lambda t: core.add(core.reduce_mod(Fi[3], p), core.scale(core.reduce_mod(Fi[5], p), t, p), p), list(range(p)), p)
    scan('F_4 + t F_5', lambda t: core.add(core.reduce_mod(Fi[4], p), core.scale(core.reduce_mod(Fi[5], p), t, p), p), list(range(p)), p)
    prs = [(t, u) for t in range(p) for u in range(1, p)]
    scan('F_3 + t F_4 + u F_5 (u != 0)', lambda tu: core.add(core.add(core.reduce_mod(Fi[3], p), core.scale(core.reduce_mod(Fi[4], p), tu[0], p), p), core.scale(core.reduce_mod(Fi[5], p), tu[1], p), p), prs, p)
# (iii) Y-even polars a F_0 + b F_1 + c F_2: non-eigen (weights off) -> unrestricted system; scan P^2(GF(101)) partially: lines c=1, b in GF, a in GF  (10201 points) 
def mk(abc, p):
    a, b, c = abc
    return core.add(core.add(core.scale(core.reduce_mod(Fi[0], p), a, p), core.scale(core.reduce_mod(Fi[1], p), b, p), p), core.scale(core.reduce_mod(Fi[2], p), c, p), p)
def scan_unres(family_name, make_D, params, p):
    t0 = time.time(); hits = []; ranks = {}
    ctx = Ctx(F, A, p)
    for t in params:
        D = make_D(t)
        ok, N0, piv = e1_lift(ctx, D, False)
        if not ok: ranks['e1fail'] = ranks.get('e1fail', 0) + 1; continue
        res = koszul_test(ctx, D, N0, False)
        key = (res['consistent'], res['rank']); ranks[str(key)] = ranks.get(str(key), 0) + 1
        if res['consistent']: hits.append(t)
    el = round(time.time() - t0, 1)
    log({'D': f'scan {family_name} over GF({p})', 'd': 5, 'family': 'scan-2param', 'engine': 'koszul', 'e1': None, 'consistent': bool(hits),
         'rank': None, 'nullity': None, 'seconds': el, 'note': f'{len(params)} parameter values; hits={hits[:20]}; rank histogram {ranks}'})
    print(f'SCAN {family_name} p={p}: hits={hits[:20]} ranks={ranks} [{el}s]', flush=True)
p = 101
pts = [(a, b, 1) for a in range(p) for b in range(p)] + [(a, 1, 0) for a in range(p)] + [(1, 0, 0)]
scan_unres('a F_0 + b F_1 + c F_2 (all of P^2(GF(101)))', lambda abc: mk(abc, p), pts, p)
