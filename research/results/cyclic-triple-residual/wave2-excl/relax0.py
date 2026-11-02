"""Degree-0 Koszul relaxation (stronger than the Hessian relaxation).
E1-admissible family D(theta) = sum_k theta_k D_k with particular E1 solutions N^(k) (N^(k)(F) = A D_k).
General E1 solution N = N_0(theta) + Koszul(s).  E2 mod F:   pi(v(theta)) + pi(W(theta)) s = 0  with
   v = D div N_0 - N_0(D)  (quadratic in theta),   W(theta) s = D div(Kz s) - (Kz s)(D)  (bilinear).
Let Wtot = span over all theta, s of the columns of pi W(theta) (= span of pi(T_S(D_k)) over basis elements), and
Y = annihilator of Wtot.  Necessary condition: y(pi v(theta)) = 0 for all y in Y  (quadrics in theta).
Modulo J(F) this is weaker (Wtot subset of pi(J)), so these quadrics contain the Hessian-relaxation quadrics.

usage: python3 relax0.py d p threads full            (d = 5..8; d = 8 needs e1_d8_exact.json)
       python3 relax0.py d p threads sector <c0>,<c1>   (d >= 9; kernel elements from kernel_exact_d{d}.json)
"""
import sys, os, time, json, re
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from w2lib import *
from multiprocessing import Pool

d = int(sys.argv[1]); p = int(sys.argv[2]); threads = int(sys.argv[3]); mode = sys.argv[4]
sector = tuple(int(x) for x in sys.argv[5].split(',')) if mode == 'sector' else None
tag = f'relax0_d{d}_p{p}' + (f'_s{sector[0]}{sector[1]}' if sector else '_full')
base = f'{PK}/{tag}'
S = Setup()
p7inv = pow(7, -1, p)   # Fz = 7 F
Fz_p = {e: v % p for e, v in S.Fz.items()}
Fz_rest = {e: v for e, v in Fz_p.items() if e[2] < 6}      # Fz - 7 X2^6

def reduce_F(P):
    """Reduce P (dict, int coeffs) modulo F over F_p: result has X2-exponent <= 5."""
    P = {e: v % p for e, v in P.items() if v % p}
    while True:
        hi = [(e, v) for e, v in P.items() if e[2] >= 6]
        if not hi: return P
        for e, v in hi:
            del P[e]
            c = (v * p7inv) % p
            m = (e[0], e[1], e[2] - 6, e[3], e[4], e[5])
            for e2, v2 in Fz_rest.items():
                k = (m[0]+e2[0], m[1]+e2[1], m[2]+e2[2], m[3]+e2[3], m[4]+e2[4], m[5]+e2[5])
                nv = (P.get(k, 0) - c * v2) % p
                if nv: P[k] = nv
                else: P.pop(k, None)

def pmul(a, b): return imul(a, b)
def padd(a, b): return core.add(a, b)
def psub(a, b): return core.sub(a, b)
def pderiv_(a, i): return ideriv(a, i)

t0 = time.time()
# ---- family: D_k and N^(k) ----
Dk = []; Nk = []; th = []
if mode == 'full':
    assert 5 <= d <= 8
    if d == 8:
        dat = json.load(open(PK + '/e1_d8_exact.json'))
        D8z = {tuple(map(int, k.split(','))): v for k, v in dat['D8z'].items()}
        Pz = [{tuple(map(int, k.split(','))): v for k, v in Pi.items()} for Pi in dat['Pz']]
        Dk.append(D8z); Nk.append(Pz); th.append(('c',))
    for i in range(6):
        for m in mons(d - 5):
            Dk.append(pmul({m: 1}, S.dF[i])); N = [{} for _ in range(6)]; N[i] = pmul(S.Az, {m: 1}); Nk.append(N); th.append((i, m))
    Sbasis = [(i, j, m) for i in range(6) for j in range(i + 1, 6) for m in mons(d - 4)]
else:
    dat = json.load(open(PK + f'/kernel_exact_d{d}.json'))
    chi_s = str(list(sector))
    for ent in dat['elements']:
        if tuple(ent['sector']) != sector: continue
        K = {tuple(map(int, k.split(','))): v for k, v in ent['Kz'].items()}
        P_ = [{tuple(map(int, k.split(','))): v for k, v in Pi.items()} for Pi in ent['Pz']]
        Dk.append(K); Nk.append(P_); th.append(('k', ent['index']))
    for i in range(6):
        for m in mons(d - 5):
            if S.chi(m) != cadd(sector, S.chiE[i]): continue
            Dk.append(pmul({m: 1}, S.dF[i])); N = [{} for _ in range(6)]; N[i] = pmul(S.Az, {m: 1}); Nk.append(N); th.append((i, m))
    target = cadd(S.chiA, sector)
    Sbasis = [(i, j, m) for i in range(6) for j in range(i + 1, 6) for m in mons(d - 4)
              if S.chi(m) == cadd(target, S.chiE[i], S.chiE[j])]
nth = len(th)
print(f'[{tag}] family size {nth}, Koszul basis {len(Sbasis)} ({time.time()-t0:.1f}s)', flush=True)

# ---- Koszul span Wtot per sector ----
dD = [[pderiv_(D, i) for i in range(6)] for D in Dk]
def koszul_vec(args):
    (i, j, m), k = args
    D = Dk[k]
    dm_i, dm_j = pderiv_({m: 1}, i), pderiv_({m: 1}, j)
    divN = psub(pmul(dm_i, S.dF[j]), pmul(dm_j, S.dF[i]))
    ND = psub(pmul(pmul({m: 1}, S.dF[j]), dD[k][i]), pmul(pmul({m: 1}, S.dF[i]), dD[k][j]))
    T = psub(pmul(D, divN), ND)
    return reduce_F(T)
def v_vec(args):
    k, l = args
    def part(k, l):
        divN = {}
        for i in range(6): divN = padd(divN, pderiv_(Nk[l][i], i))
        ND = {}
        for i in range(6): ND = padd(ND, pmul(Nk[l][i], dD[k][i]))
        return psub(pmul(Dk[k], divN), ND)
    T = part(k, l)
    if k != l: T = padd(T, part(l, k))
    return reduce_F(T)
tasks_w = [((i, j, m), k) for (i, j, m) in Sbasis for k in range(nth)]
tasks_v = [(k, l) for k in range(nth) for l in range(k, nth)]
t1 = time.time()
with Pool(threads) as pool:
    Wv = pool.map(koszul_vec, tasks_w, chunksize=64)
    Vv = pool.map(v_vec, tasks_v, chunksize=64)
print(f'[{tag}] computed {len(Wv)} Koszul images and {len(Vv)} v-pairs ({time.time()-t1:.1f}s)', flush=True)
# group by character (character of the polynomial; zero polynomials skipped)
bych = {}
for w in Wv:
    if w: bych.setdefault(chi_of(w, S.chi), [[], []])[0].append(w)
for (k, l), v in zip(tasks_v, Vv):
    if v: bych.setdefault(chi_of(v, S.chi), [[], []])[1].append(((k, l), v))
quadrics = []; info = {}
for c, (ws, vs) in sorted(bych.items()):
    if not vs: continue
    monsc = sorted({e for w in ws for e in w} | {e for _, v in vs for e in v})
    idx = {e: i for i, e in enumerate(monsc)}; C = len(monsc)
    # Wtot_c: rows = Koszul images; annihilator Y = left-null... we need functionals y on the monomial space with y(w) = 0 for all w:
    # i.e. nullspace of the matrix with rows = w vectors.
    if ws:
        Mw = nmod_mat(len(ws), C, [0] * (len(ws) * C), p)
        flat = [0] * (len(ws) * C)
        for r, w in enumerate(ws):
            for e, v in w.items(): flat[r * C + idx[e]] = v % p
        Mw = nmod_mat(len(ws), C, flat, p)
        Yn, ny = Mw.nullspace()          # columns of Yn (C x ny) are the functionals
        rkW = C - ny if ny is not None else C
        Y = nmod_mat(C, ny, [int(Yn[i, j]) for i in range(C) for j in range(ny)], p) if ny else None
    else:
        rkW = 0; ny = C
        Y = nmod_mat(C, C, [int(i == j) for i in range(C) for j in range(C)], p)
    if ny == 0:
        info[str(c)] = {'mons': C, 'rank_W': rkW, 'n_functionals': 0, 'n_vpairs': len(vs)}; continue
    flatv = [0] * (len(vs) * C)
    for r, (_, v) in enumerate(vs):
        for e, val in v.items(): flatv[r * C + idx[e]] = val % p
    Mv = nmod_mat(len(vs), C, flatv, p)
    Q = Mv * Y                              # (n_vpairs x ny): coefficient of pair r in functional j
    nq = 0
    for j in range(ny):
        qd = {}
        for r, ((k, l), _) in enumerate(vs):
            val = int(Q[r, j])
            if val: qd[(k, l)] = val
        if qd: quadrics.append(qd); nq += 1
    info[str(c)] = {'mons': C, 'rank_W': rkW, 'n_functionals': ny, 'n_vpairs': len(vs), 'n_quadrics': nq}
print(f'[{tag}] quadrics: {len(quadrics)}; sectors: {info} ({time.time()-t0:.1f}s)', flush=True)
names = [f't{k}' for k in range(nth)]
polys = [quad_str(q, names, p) for q in quadrics]
write_msolve(base + '.ms', nth, p, polys, names)
json.dump({'d': d, 'p': p, 'mode': mode, 'sector': sector, 'theta': [list(t) if t[0] not in ('c', 'k') else list(map(str, t)) for t in th],
           'n_koszul': len(Sbasis), 'n_quadrics': len(quadrics), 'sectors': info}, open(base + '_system.json', 'w'))
log({'exp': 'koszul-relaxation-build', 'family': f'E1-admissible degree {d} ' + (f'sector {sector}' if sector else 'all'), 'd': d, 'prime': p,
     'result': f'{len(quadrics)} quadrics in {nth} unknowns (Koszul span ranks per sector in receipt)', 'rank': None,
     'seconds': round(time.time() - t0, 1), 'job': f'excl-ctr-{tag}', 'note': 'degree-0 Koszul relaxation: functionals killing F*S and all Koszul images; modular'})
# ---- Groebner basis + illegal-locus forms ----
rc, dt, so, se = run_msolve(base + '.ms', base + '_gb.txt', threads=threads, gb=2)
gb = parse_msolve_gb(base + '_gb.txt') if rc == 0 else []
def is_linear(s):
    s2 = s.replace(' ', '')
    return not re.search(r't\d+\^', s2) and not re.search(r't\d+\*t\d+', s2)
nlin = sum(1 for g in gb if is_linear(g))
print(f'[{tag}] GB: rc={rc} {len(gb)} elements, {nlin} linear ({dt:.1f}s)', flush=True)
log({'exp': 'koszul-relaxation-gb', 'family': f'E1-admissible degree {d} ' + (f'sector {sector}' if sector else 'all'), 'd': d, 'prime': p,
     'result': f'GB has {len(gb)} elements, {nlin} linear', 'rank': None, 'seconds': round(dt, 1), 'job': f'excl-ctr-{tag}', 'note': f'rc={rc}; msolve -g 2 mod p'})
# gauge-fixed m-primary test (rigorous; see relax_run.py): slice theta_(5, m) = 0 for Y2 | m
removed = {k for k, t in enumerate(th) if t[0] == 5 and t[1][5] >= 1}
keep = [k for k in range(nth) if k not in removed]
kidx = {k: i for i, k in enumerate(keep)}
gq = []
for q in quadrics:
    q2 = {(kidx[a], kidx[b]): v for (a, b), v in q.items() if a not in removed and b not in removed}
    if q2: gq.append(q2)
gnames = [f'u{i}' for i in range(len(keep))]
gpolys = [quad_str(q, gnames, p) for q in gq]
write_msolve(base + '_gauge.ms', len(keep), p, gpolys, gnames)
rc_g, dt_g, so_g, se_g = run_msolve(base + '_gauge.ms', base + '_gauge_lead.txt', threads=threads, gb=1)
lead = parse_msolve_gb(base + '_gauge_lead.txt') if rc_g == 0 else []
pure = set()
for mon in lead:
    vs = set(re.findall(r'u\d+', mon))
    if len(vs) == 1: pure.add(vs.pop())
m_primary = (rc_g == 0) and all(f'u{i}' in pure for i in range(len(keep)))
print(f'[{tag}] gauge-fixed: {len(keep)} vars ({len(removed)} removed), {len(gq)} quadrics; pure powers {len(pure)}/{len(keep)} => m-primary={m_primary} ({dt_g:.1f}s)', flush=True)
log({'exp': 'koszul-relaxation-gauge-fixed-mprimary', 'family': f'E1-admissible degree {d} ' + (f'sector {sector}' if sector else 'all'), 'd': d, 'prime': p,
     'result': ('m-PRIMARY mod p => no denominator in this family (rigorous via properness)' if m_primary else f'NOT m-primary (pure powers {len(pure)}/{len(keep)})'),
     'rank': None, 'seconds': round(dt_g, 1), 'job': f'excl-ctr-{tag}', 'note': f'{len(gq)} quadrics in {len(keep)} unknowns after gauge; msolve -g 1; rc={rc_g}'})
json.dump({'n_vars_gauge': len(keep), 'n_removed': len(removed), 'n_quadrics': len(gq), 'lead_size': len(lead), 'pure_powers': len(pure), 'm_primary': m_primary,
           'seconds': round(dt_g, 1), 'rc': rc_g}, open(base + '_gauge.json', 'w'), indent=1)
if os.environ.get('SKIP_FORMS'):
    sys.exit(0)
# illegal forms
idx = {t: k for k, t in enumerate(th)}
forms = []
for t in th:
    if t[0] in ('c', 'k'): forms.append({idx[t]: 1})
    elif t[1][t[0]] == 0: forms.append({idx[t]: 1})
if d >= 6:
    for m in mons(d - 6):
        if sector is not None and S.chi(m) != sector: continue
        es = []
        for i in range(6):
            e = list(m); e[i] += 1; e = tuple(e)
            if (i, e) in idx: es.append(idx[(i, e)])
        for j in range(1, len(es)):
            forms.append({es[0]: 1, es[j]: -1})
os.makedirs(base + '_forms', exist_ok=True)
def one(k):
    f = forms[k]
    inp = f'{base}_forms/f{k}.ms'; out = f'{base}_forms/f{k}.out'
    write_msolve(inp, nth, p, polys + [lin_str(f, names) + '-1'], names)
    rc, dt, so, se = run_msolve(inp, out, threads=1, gb=0)
    txt = open(out).read().strip() if os.path.exists(out) else ''
    return {'k': k, 'form': lin_str(f, names), 'rc': rc, 'seconds': round(dt, 1), 'empty': txt.startswith('[-1]'), 'out': txt[:200]}
t1 = time.time()
with Pool(min(len(forms), threads)) as pool:
    res = pool.map(one, range(len(forms)))
n_empty = sum(1 for r in res if r['empty'])
json.dump({'d': d, 'p': p, 'mode': mode, 'sector': sector, 'n_theta': nth, 'n_quadrics': len(quadrics), 'gb_size': len(gb), 'gb_linear': nlin,
           'n_forms': len(forms), 'n_empty': n_empty, 'forms': res, 'seconds_forms': round(time.time() - t1, 1), 'sectors': info},
          open(base + '.json', 'w'), indent=1)
print(f'[{tag}] forms: {n_empty}/{len(forms)} empty ({time.time()-t1:.1f}s)', flush=True)
log({'exp': 'koszul-relaxation-illegal-locus', 'family': f'E1-admissible degree {d} ' + (f'sector {sector}' if sector else 'all'), 'd': d, 'prime': p,
     'result': f'{n_empty}/{len(forms)} illegal-locus forms l have V(I+(l-1)) empty mod p' + (' => V(I) inside illegal locus (mod p)' if n_empty == len(forms) else ' => EXTRA components'),
     'rank': None, 'seconds': round(time.time() - t1, 1), 'job': f'excl-ctr-{tag}', 'note': 'msolve emptiness per form; modular evidence'})
