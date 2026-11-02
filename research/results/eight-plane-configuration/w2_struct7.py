"""Evaluate the exact d=7 Hessian quadrics (mod p) at structured quadratic vector fields."""
import sys, json, itertools, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, mul, add, scale, pderiv, parse
from fractions import Fraction
WD = '/tmp/hodge/results/eight-plane-configuration'
P = core.PRIMES[1]
def lin(i):
    e = [0]*6; e[i] = 1; return {tuple(e): Fraction(1)}
Z = [lin(i) for i in range(6)]
Q = [add(mul(Z[i], Z[i]), mul(Z[i+3], Z[i+3])) for i in range(3)]
# linear fields as lists of 6 polys
def field(comps): return comps
E = [Z[i] for i in range(6)]
EX = [Z[0], Z[1], Z[2], {}, {}, {}]; EY = [{}, {}, {}, Z[3], Z[4], Z[5]]
S = [[Z[j] if j in (i, i+3) else {} for j in range(6)] for i in range(3)]
R = [[(scale(Z[i+3], -1) if j == i else (Z[i] if j == i+3 else {})) for j in range(6)] for i in range(3)]
basic = {'E': E, 'EX': EX, 'EY': EY}
for i in range(3): basic[f'S{i}'] = S[i]; basic[f'R{i}'] = R[i]
for a in range(6):
    for b in range(6):
        basic[f'z{a}d{b}'] = [Z[a] if j == b else {} for j in range(6)]
cands = {}
# linear form * linear field
for a in range(6):
    for name, v in basic.items():
        cands[f'z{a}*{name}'] = [mul(Z[a], c) if c else {} for c in v]
# Q_j * constant field
for j in range(3):
    for b in range(6):
        cands[f'Q{j}*d{b}'] = [Q[j] if k == b else {} for k in range(6)]
# gradient fields of cubics: monomials, z_a Q_j, and sums
cubics = {}
for m in monomials(3, None):
    cubics['grad(' + core.to_str({m: Fraction(1)}) + ')'] = {m: Fraction(1)}
for a in range(6):
    for j in range(3):
        cubics[f'grad(z{a}Q{j})'] = mul(Z[a], Q[j])
for name, G in cubics.items():
    cands[name] = [pderiv(G, i) for i in range(6)]
# Hamiltonian-type: F_i-free quadratic fields built from Q's: Q_j d_a - ... skip
results = []
for par in (0, 1):
    dat = json.load(open(f'{WD}/hess7local-quadrics-d7-p{par}.json'))
    Wm = [(i, tuple(m)) for i, m in dat['Wm']]
    widx = {wm: k for k, wm in enumerate(Wm)}
    rows = [{(k, l): v for (k, l), v in r} for r in dat['rows']]
    for name, comps in cands.items():
        # coordinate vector (mod p) in the Wm basis; require pure parity = par
        coords = {}
        ok = True
        for i in range(6):
            for m, c in comps[i].items():
                if (i, m) not in widx:
                    ok = False; break
                cc = Fraction(c); coords[widx[(i, m)]] = (cc.numerator * pow(cc.denominator, -1, P)) % P
            if not ok: break
        if not ok or not coords:
            continue
        nz = 0
        for row in rows:
            s = 0
            for (k, l), v in row.items():
                if k in coords and l in coords:
                    s += v * coords[k] * coords[l]
            if s % P: nz += 1
        results.append({'d': 7, 'par': par, 'field': name, 'nonzero_quadrics': nz, 'of': len(rows)})
hits = [r for r in results if r['nonzero_quadrics'] == 0]
print('evaluated', len(results), 'structured fields; Hessian-admissible:', len(hits))
for r in hits: print(r)
zs = sorted(results, key=lambda r: r['nonzero_quadrics'])[:8]
print('closest (fewest nonzero quadrics):', zs)
json.dump(results, open(f'{WD}/struct7-results.json', 'w'))
with open(f'{WD}/wave2-search.log.jsonl', 'a') as f:
    f.write(json.dumps({'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'structured quadratic fields vs d=7 Hessian quadrics (local, mod p)', 'd': 7, 'n_fields': len(results), 'admissible': [r['field'] for r in hits], 'closest': zs[:5]}) + '\n')
