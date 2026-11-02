#!/usr/bin/env python3
"""w2_points.py -- exact algebraic solving of the stage-1 system at a Jacobian-type degree with msolve.

Builds the quadratic forms q_r(v) (v = quotient w-coordinates, see w2_closure.py) over F_p (--p) or over Q (--char0,
multi-modular normal forms + rational reconstruction), writes an msolve file with a random dehomogenizing linear
equation l(v) = 1, runs msolve, and reports the solution structure (dimension / degree / real roots).
usage: python3 w2_points.py --d 6 --parity 0 [--p 2147483629] [--char0] [--threads 8] [--timeout 3600]
"""
import sys, os, json, time, argparse, itertools, subprocess, random
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from fractions import Fraction
from flint import nmod_mat
from hodge import core, jacobian
from hodge.core import NV, monomials, mul, add, scale, pderiv, reduce_mod, shift, sub

LANE = 'quadratic-cubic-pullback'
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True); ap.add_argument('--parity', type=int, default=0)
ap.add_argument('--p', type=int, default=2147483629); ap.add_argument('--char0', action='store_true')
ap.add_argument('--nprimes', type=int, default=12, help='primes for the char-0 reconstruction')
ap.add_argument('--threads', type=int, default=4); ap.add_argument('--timeout', type=int, default=3600)
ap.add_argument('--instance', default=LANE); ap.add_argument('--fermat', action='store_true')
ap.add_argument('--out', default=None); ap.add_argument('--no-dehom', action='store_true')
a = ap.parse_args()
d, parD = a.d, a.parity
F, A, q, _ = core.load_instance(a.instance)
if a.fermat:
    F = {tuple(6 if j == i else 0 for j in range(6)): Fraction(1) for i in range(6)}
t0 = time.time()

def forms_mod(p):
    Fp, Ap = reduce_mod(F, p), reduce_mod(A, p)
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    wlab = []
    for l in range(NV):
        par_w = (parD - (0 if l < 3 else 1)) % 2
        for mu in monomials(d - 5, par_w):
            wlab.append((l, mu))
    n = len(wlab)
    Dcols = [shift(Fi[l], mu, 1, p) for (l, mu) in wlab]
    Ncols = [(l, shift(Ap, mu, 1, p)) for (l, mu) in wlab]
    widx = {lab: i for i, lab in enumerate(wlab)}
    spur = []
    inv6 = pow(6, -1, p)
    for G in (monomials(d - 6, parD) if d >= 6 else []):
        v = [0] * n
        for l in range(NV):
            mu = list(G); mu[l] += 1; v[widx[(l, tuple(mu))]] = inv6
        spur.append(v)
    m = len(spur)
    pivots = []
    if m:
        Sm = nmod_mat(m, n, [x for row in spur for x in row], p)
        Sr, rk = Sm.rref(); assert rk == m
        for i in range(rk):
            j = 0
            while int(Sr[i, j]) == 0: j += 1
            pivots.append(j)
    free = [j for j in range(n) if j not in set(pivots)]
    red = jacobian._Reducer(Fp, 2 * d, p, 1)
    sm = red.standard; sidx = {mo: i for i, mo in enumerate(sm)}
    def Bpair(a_, b_):
        Da = Dcols[a_]; lb, Nb = Ncols[b_]
        out = mul(Da, pderiv(Nb, lb, p), p)
        return sub(out, mul(Nb, pderiv(Da, lb, p), p), p)
    nq = len(free)
    pairs = [(x, y) for x in range(nq) for y in range(x, nq)]
    rows = [[0] * len(pairs) for _ in range(len(sm))]
    for j, (x, y) in enumerate(pairs):
        P = Bpair(free[x], free[y])
        if x != y:
            P = add(P, Bpair(free[y], free[x]), p)
        for mo, c in red.normal_form(P).items():
            rows[sidx[mo]][j] = c
    M = nmod_mat(len(sm), len(pairs), [x for row in rows for x in row], p)
    Mr, rk = M.rref()
    forms = []
    for i in range(rk):
        forms.append({pairs[j]: int(Mr[i, j]) for j in range(len(pairs)) if int(Mr[i, j])})
    return forms, nq, pairs, wlab, free

if not a.char0:
    p = a.p
    forms, nq, pairs, wlab, free = forms_mod(p)
    char = p
    def cstr(c):
        return str(c if c <= p // 2 else c - p)
else:
    # multi-modular reconstruction of the reduced forms over Q (rref pivot structure assumed stable across primes)
    primes = core.PRIMES[:a.nprimes]
    acc = None
    for k, p in enumerate(primes):
        forms, nq, pairs, wlab, free = forms_mod(p)
        vec = {}
        for i, fr in enumerate(forms):
            for key, c in fr.items():
                vec[(i, key)] = c
        if acc is None:
            acc = ({kk: v for kk, v in vec.items()}, p)
        else:
            xs, mmod = acc
            keys = set(xs) | set(vec)
            xs = {kk: core.crt_pair(xs.get(kk, 0), mmod, vec.get(kk, 0), p)[0] for kk in keys}
            acc = (xs, mmod * p)
        xs, mmod = acc
        rat = {kk: core.ratrecon(v, mmod) for kk, v in xs.items()}
        if all(v is not None for v in rat.values()):
            print(f'[points] rational reconstruction stable after {k+1} primes', flush=True)
            break
    forms = [{} for _ in range(len(forms))]
    for (i, key), v in rat.items():
        if v is None:
            raise SystemExit('rational reconstruction failed; increase --nprimes')
        if v != 0:
            forms[i][key] = v
    char = 0
    def cstr(c):
        return str(Fraction(c))

print(f'[points d={d} par={parD} char={char}] {len(forms)} quadratic forms in {nq} variables ({time.time()-t0:.1f}s)', flush=True)
names = [f'v{i}' for i in range(nq)]
lines = [','.join(names), str(char)]
polys = []
for fr in forms:
    terms = []
    for (x, y), c in fr.items():
        mon = f'{names[x]}^2' if x == y else f'{names[x]}*{names[y]}'
        terms.append(f'{cstr(c)}*{mon}')
    polys.append('+'.join(terms).replace('+-', '-'))
if not a.no_dehom:
    rng = random.Random(1234)
    lin = '+'.join(f'{rng.randint(1, 97)}*{nm}' for nm in names) + '-1'
    polys.append(lin)
out = a.out or f'/tmp/hodge/results/{LANE}/w2_points_d{d}_p{parD}_c{char}.ms'
with open(out, 'w') as f:
    f.write('\n'.join(lines) + '\n' + ',\n'.join(polys) + '\n')
print(f'[points] wrote {out}', flush=True)
res_path = out.replace('.ms', '.out')
cmd = ['msolve', '-t', str(a.threads), '-v', '1', '-f', out, '-o', res_path]
if char == 0:
    cmd += ['-P', '1']
t1 = time.time()
try:
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
    print(f'[points] msolve exit {r.returncode} in {time.time()-t1:.1f}s', flush=True)
    print(r.stderr[-3000:], flush=True)
    print(r.stdout[-2000:], flush=True)
except subprocess.TimeoutExpired:
    print(f'[points] msolve TIMEOUT after {a.timeout}s', flush=True)
if os.path.exists(res_path):
    txt = open(res_path).read()
    print(f'[points] output {res_path}: {len(txt)} chars; head: {txt[:1500]}', flush=True)
rec = {'lane': LANE, 'item': 'points', 'd': d, 'parity': parD, 'char': char, 'n_forms': len(forms), 'nq': nq,
       'msolve_seconds': round(time.time() - t1, 1), 'out': res_path, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
if os.path.exists(res_path):
    rec['out_head'] = open(res_path).read()[:400]
with open(f'/tmp/hodge/results/{LANE}/w2_points.jsonl', 'a') as f:
    f.write(json.dumps(rec) + '\n')
