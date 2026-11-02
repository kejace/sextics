"""Full parametric E2 system for the E1-admissible family of a character sector at degree d (semi-invariant D):
unknowns theta (family coefficients) and s (Koszul coefficients); equations: coefficients of
   pi( D(theta) div N_0(theta) - N_0(theta)(D(theta)) + D div(Kz s) - (Kz s)(D) )  = 0,   pi = reduction mod F,
a system of quadrics (quadratic in theta, bilinear in (theta, s)).  Its solutions with theta outside the illegal locus
L = {c = 0, G = h E} are exactly the degree-d certificates with denominator in the sector.  We compute the reduced
Groebner basis mod p and run per-form emptiness tests  I + (l(theta) - 1)  for the linear forms l cutting out L
(MODULAR EVIDENCE: the system is not gauge-invariant, so no properness argument applies).
usage: python3 bilin.py d p threads c0,c1 [timeout_s]"""
import sys, os, time, json, re
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from multiprocessing import Pool
d = int(sys.argv[1]); p = int(sys.argv[2]); threads = int(sys.argv[3]); sector = tuple(int(x) for x in sys.argv[4].split(','))
tmo = int(sys.argv[5]) if len(sys.argv) > 5 else 7200
tag = f'bilin_d{d}_p{p}_s{sector[0]}{sector[1]}'; base = f'{PK}/{tag}'
S = Setup()
p7inv = pow(7, -1, p); Fz_rest = {e: v % p for e, v in S.Fz.items() if e[2] < 6}
def reduce_F(P):
    P = {e: v % p for e, v in P.items() if v % p}
    while True:
        hi = sorted((e for e in P if e[2] >= 6), key=lambda e: -e[2])
        if not hi: return P
        for e in hi:
            v = P.pop(e, None)
            if v is None: continue
            c = (v * p7inv) % p; m = (e[0], e[1], e[2] - 6, e[3], e[4], e[5])
            for e2, v2 in Fz_rest.items():
                k = (m[0]+e2[0], m[1]+e2[1], m[2]+e2[2], m[3]+e2[3], m[4]+e2[4], m[5]+e2[5])
                nv = (P.get(k, 0) - c * v2) % p
                if nv: P[k] = nv
                else: P.pop(k, None)
t0 = time.time()
Dk = []; Nk = []; th = []
# kernel elements: exact file if available, else mod p canonical
kf = PK + f'/kernel_exact_d{d}.json'
nker = 0
if os.path.exists(kf):
    dat = json.load(open(kf))
    for ent in dat['elements']:
        if tuple(ent['sector']) != sector or not ent.get('verified'): continue
        Dk.append({tuple(map(int, k.split(','))): v for k, v in ent['Kz'].items()})
        Nk.append([{tuple(map(int, k.split(','))): v for k, v in Pi.items()} for Pi in ent['Pz']]); th.append(('k', ent['index'])); nker += 1
if nker == 0:
    Dm, vecs, rJ = kernel_canonical(S, d, sector, p)
    target = cadd(S.chiA, sector)
    mT = [m for m in mons(d + 6) if S.chi(m) == target]; iT = {m: i for i, m in enumerate(mT)}
    colsE = []; metaE = []
    for i in range(6):
        ci = cadd(target, S.chiE[i])
        for m in mons(d + 1):
            if S.chi(m) != ci: continue
            colsE.append(imul({m: 1}, S.dF[i])); metaE.append((i, m))
    for v in vecs:
        K = {Dm[j]: v[j] for j in range(len(Dm)) if v[j]}
        rhs = imul(S.Az, K)
        R, C = len(mT), len(colsE)
        flat = [0] * (R * (C + 1))
        for j, P_ in enumerate(colsE):
            for e, val in P_.items(): flat[iT[e] * (C + 1) + j] = (flat[iT[e] * (C + 1) + j] + val) % p
        for e, val in rhs.items(): flat[iT[e] * (C + 1) + C] = val % p
        Mr, rk = nmod_mat(R, C + 1, flat, p).rref()
        piv = []
        for i in range(rk):
            j = 0
            while j <= C and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        assert piv[-1] != C, 'E1 inconsistent for a kernel element?!'
        N = [{} for _ in range(6)]
        for i, j in enumerate(piv):
            val = int(Mr[i, C])
            if val: N[metaE[j][0]][metaE[j][1]] = val
        Dk.append(K); Nk.append(N); th.append(('k', len(th))); nker += 1
    print(f'kernel elements mod p: {nker}', flush=True)
for i in range(6):
    for m in mons(d - 5):
        if S.chi(m) != cadd(sector, S.chiE[i]): continue
        Dk.append(imul({m: 1}, S.dF[i])); N = [{} for _ in range(6)]; N[i] = imul(S.Az, {m: 1}); Nk.append(N); th.append((i, m))
nth = len(th)
target = cadd(S.chiA, sector)
Sbasis = [(i, j, m) for i in range(6) for j in range(i + 1, 6) for m in mons(d - 4) if S.chi(m) == cadd(target, S.chiE[i], S.chiE[j])]
nS = len(Sbasis)
print(f'[{tag}] theta {nth} (kernel {nker}), Koszul {nS} ({time.time()-t0:.1f}s)', flush=True)
dD = [[ideriv(D, i) for i in range(6)] for D in Dk]
def kvec(args):
    (i, j, m), k = args
    D = Dk[k]
    divN = core.sub(imul(ideriv({m: 1}, i), S.dF[j]), imul(ideriv({m: 1}, j), S.dF[i]))
    ND = core.sub(imul(imul({m: 1}, S.dF[j]), dD[k][i]), imul(imul({m: 1}, S.dF[i]), dD[k][j]))
    return reduce_F(core.sub(imul(D, divN), ND))
def vvec(args):
    k, l = args
    def part(k, l):
        divN = {}
        for i in range(6): divN = core.add(divN, ideriv(Nk[l][i], i))
        ND = {}
        for i in range(6): ND = core.add(ND, imul(Nk[l][i], dD[k][i]))
        return core.sub(imul(Dk[k], divN), ND)
    T = part(k, l)
    if k != l: T = core.add(T, part(l, k))
    return reduce_F(T)
tw = [(s_, k) for s_ in Sbasis for k in range(nth)]; tv = [(k, l) for k in range(nth) for l in range(k, nth)]
with Pool(threads) as pool:
    Wv = pool.map(kvec, tw, chunksize=32); Vv = pool.map(vvec, tv, chunksize=32)
# equations indexed by monomials
eqs = {}
for (k, l), v in zip(tv, Vv):
    for e, val in v.items(): eqs.setdefault(e, {})[('t', k, l)] = val
for (s_, k), w in zip(tw, Wv):
    si = Sbasis.index(s_)
    for e, val in w.items(): eqs.setdefault(e, {})[('s', si, k)] = val
names = [f't{k}' for k in range(nth)] + [f's{i}' for i in range(nS)]
polys = []
for e, terms in eqs.items():
    parts = []
    for key, val in terms.items():
        val %= p
        if val > p // 2: val -= p
        if val == 0: continue
        if key[0] == 't':
            mon = f't{key[1]}^2' if key[1] == key[2] else f't{key[1]}*t{key[2]}'
        else:
            mon = f's{key[1]}*t{key[2]}'
        parts.append(f'{"+" if val > 0 else "-"}{abs(val)}*{mon}')
    if parts: polys.append(''.join(parts).lstrip('+'))
print(f'[{tag}] {len(polys)} equations in {nth + nS} unknowns ({time.time()-t0:.1f}s)', flush=True)
write_msolve(base + '.ms', nth + nS, p, polys, names)
log({'exp': 'bilinear-sector-system-build', 'family': f'degree {d} sector {sector}', 'd': d, 'prime': p, 'result': f'{len(polys)} quadratic equations in {nth} theta + {nS} s unknowns',
     'rank': None, 'seconds': round(time.time() - t0, 1), 'job': f'excl-ctr-{tag}', 'note': 'full parametric E2 for semi-invariant D (E1 pre-solved); modular'})
rc, dt, so, se = run_msolve(base + '.ms', base + '_gb.txt', threads=threads, gb=2, timeout=tmo) if True else (1, 0, '', '')
gb = parse_msolve_gb(base + '_gb.txt') if rc == 0 else []
print(f'[{tag}] GB rc={rc} size {len(gb)} ({dt:.1f}s)', flush=True)
log({'exp': 'bilinear-sector-system-gb', 'family': f'degree {d} sector {sector}', 'd': d, 'prime': p, 'result': f'reduced GB size {len(gb)} (rc={rc})', 'rank': None,
     'seconds': round(dt, 1), 'job': f'excl-ctr-{tag}', 'note': 'msolve -g 2 mod p'})
# per-form emptiness tests on theta-forms cutting out the illegal locus
idx = {t: k for k, t in enumerate(th)}
forms = [{idx[t]: 1} for t in th if t[0] == 'k' or t[1][t[0]] == 0]
for m in mons(d - 6):
    if S.chi(m) != sector: continue
    es = []
    for i in range(6):
        e = list(m); e[i] += 1; e = tuple(e)
        if (i, e) in idx: es.append(idx[(i, e)])
    for j in range(1, len(es)): forms.append({es[0]: 1, es[j]: -1})
os.makedirs(base + '_forms', exist_ok=True)
src_polys = gb if (rc == 0 and gb) else polys
def one(k):
    f = forms[k]; inp = f'{base}_forms/f{k}.ms'; out = f'{base}_forms/f{k}.out'
    write_msolve(inp, nth + nS, p, src_polys + [lin_str(f, names) + '-1'], names)
    rc_, dt_, so_, se_ = run_msolve(inp, out, threads=1, gb=0, timeout=tmo)
    txt = open(out).read().strip() if os.path.exists(out) else ''
    return {'k': k, 'form': lin_str(f, names), 'rc': rc_, 'seconds': round(dt_, 1), 'empty': txt.startswith('[-1]'), 'out': txt[:120]}
t1 = time.time()
with Pool(min(len(forms), threads)) as pool:
    res = pool.map(one, range(len(forms)))
n_empty = sum(1 for r in res if r['empty'])
json.dump({'d': d, 'p': p, 'sector': sector, 'n_theta': nth, 'n_kernel': nker, 'n_s': nS, 'n_eqs': len(polys), 'gb_size': len(gb), 'gb_rc': rc,
           'n_forms': len(forms), 'n_empty': n_empty, 'forms': res}, open(base + '.json', 'w'), indent=1)
print(f'[{tag}] forms: {n_empty}/{len(forms)} empty ({time.time()-t1:.1f}s)', flush=True)
log({'exp': 'bilinear-sector-system-forms', 'family': f'degree {d} sector {sector}', 'd': d, 'prime': p,
     'result': f'{n_empty}/{len(forms)} illegal-locus forms give empty V(I+(l-1)) mod p' + (' => no semi-invariant certificate in this sector (MODULAR EVIDENCE)' if n_empty == len(forms) else ' => candidate component!'),
     'rank': None, 'seconds': round(time.time() - t1, 1), 'job': f'excl-ctr-{tag}', 'note': 'per-form msolve emptiness; not gauge-invariant, hence modular evidence only'})
