"""Evaluate exact Hessian quadrics (mod p, from <tag>-quadrics-d{d}-p{par}.json) at structured vector fields of
component degree d-5.  Candidates: (monomial of degree d-6) x basic linear fields (Euler, X/Y-Euler, block scalings,
block rotations, coordinate fields z_a d_b), Q_j x (monomial of degree d-7) x basic fields, gradients of monomials of
degree d-4 and of z^m * Q_j (deg m = d-6), and Q_j Q_k x constant fields when d = 9.

usage: python3 w2_structprobe.py --d 8 --tag hessB [--pars 0,1]
"""
import sys, json, itertools, time, argparse, os
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--tag', default='hessB')
ap.add_argument('--pars', default='0,1')
args = ap.parse_args()
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, mul, add, scale, pderiv
from fractions import Fraction
WD = '/tmp/hodge/results/eight-plane-configuration'
P = core.PRIMES[1]
d = args.d; k = d - 5
def lin(i):
    e = [0] * 6; e[i] = 1; return {tuple(e): Fraction(1)}
Z = [lin(i) for i in range(6)]
Q = [add(mul(Z[i], Z[i]), mul(Z[i + 3], Z[i + 3])) for i in range(3)]
def mono(m): return {tuple(m): Fraction(1)}
E = [Z[i] for i in range(6)]
EX = [Z[0], Z[1], Z[2], {}, {}, {}]; EY = [{}, {}, {}, Z[3], Z[4], Z[5]]
S = [[Z[j] if j in (i, i+3) else {} for j in range(6)] for i in range(3)]
R = [[(scale(Z[i + 3], -1) if j == i else (Z[i] if j == i + 3 else {})) for j in range(6)] for i in range(3)]
basic = {'E': E, 'EX': EX, 'EY': EY}
for i in range(3): basic[f'S{i}'] = S[i]; basic[f'R{i}'] = R[i]
for a in range(6):
    for b in range(6):
        basic[f'z{a}d{b}'] = [Z[a] if j == b else {} for j in range(6)]
cands = {}
for m in monomials(k - 1, None):
    for name, v in basic.items():
        cands[f'{core.to_str(mono(m))}*{name}'] = [mul(mono(m), c) if c else {} for c in v]
if k >= 2:
    for j in range(3):
        for m in monomials(k - 2, None):
            for name, v in basic.items():
                if name.startswith('z') or name in ('E', 'EX', 'EY') or name.startswith('S') or name.startswith('R'):
                    cands[f'Q{j}*{core.to_str(mono(m))}*{name}'] = [mul(mul(Q[j], mono(m)), c) if c else {} for c in v]
for m in monomials(k + 1, None):
    G = mono(m); cands['grad(' + core.to_str(G) + ')'] = [pderiv(G, i) for i in range(6)]
for j in range(3):
    for m in monomials(k - 1, None):
        G = mul(mono(m), Q[j]); cands[f'grad({core.to_str(mono(m))}*Q{j})'] = [pderiv(G, i) for i in range(6)]
if k >= 4:
    for j in range(3):
        for l in range(j, 3):
            for b in range(6):
                cands[f'Q{j}Q{l}*d{b}'] = [mul(Q[j], Q[l]) if i == b else {} for i in range(6)]
print('candidates:', len(cands), flush=True)
results = []
for par in [int(x) for x in args.pars.split(',')]:
    path = f'{WD}/{args.tag}-quadrics-d{d}-p{par}.json'
    if not os.path.exists(path):
        print('missing', path); continue
    dat = json.load(open(path))
    Wm = [(i, tuple(m)) for i, m in dat['Wm']]; widx = {wm: kk for kk, wm in enumerate(Wm)}
    rows = [{(a, b): v for (a, b), v in r} for r in dat['rows']]
    for name, comps in cands.items():
        coords = {}; ok = True
        for i in range(6):
            for m, c in comps[i].items():
                if (i, m) not in widx: ok = False; break
                cc = Fraction(c); coords[widx[(i, m)]] = (cc.numerator * pow(cc.denominator, -1, P)) % P
            if not ok: break
        if not ok or not coords: continue
        nz = 0
        for row in rows:
            s = 0
            for (a, b), v in row.items():
                if a in coords and b in coords: s += v * coords[a] * coords[b]
            if s % P: nz += 1
        results.append({'d': d, 'par': par, 'field': name, 'nonzero_quadrics': nz, 'of': len(rows)})
hits = [r for r in results if r['nonzero_quadrics'] == 0]
print('evaluated', len(results), 'fields; Hessian-admissible:', len(hits))
for r in hits: print(r)
json.dump(results, open(f'{WD}/structprobe-d{d}-results.json', 'w'))
with open(f'{WD}/wave2-search.log.jsonl', 'a') as f:
    f.write(json.dumps({'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': f'structured fields vs d={d} Hessian quadrics (mod p)', 'd': d, 'n_fields': len(results),
                        'admissible': [r['field'] for r in hits][:50], 'n_admissible': len(hits)}) + '\n')
