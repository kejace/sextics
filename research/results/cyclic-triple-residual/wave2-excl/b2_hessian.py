"""(b) Hessian relaxation of E2 for E1-admissible denominators, as a quadratic system, with m-primary test.

Setting (q = 1).  Every E1-admissible D of degree d is  D = sum_k lam_k Dk + sum_i w_i F_i  with Dk a basis of
Ann([A])_d mod J(F) (exact E1 lifts Nk: sum_i Nk_i F_i = A Dk) and w_i in S_{d-5}.  N = sum_k lam_k Nk + A w + K
(K Koszul) is the general E1 lift.  E2 requires  B_D(N) = D div N - N(D) in F*S + im(T_D) subset J(F).
Hence  Phi(lam, w) := [B_D(sum lam_k Nk + A w)] in R_{2d} = S/J  must vanish: a system of quadrics in (lam, w)
whose solution set contains the Euler gauge  w -> w + l Z/6  (D -> D + l F).  Fixing the gauge by dropping the
coordinates (Y2-component, monomials divisible by Y2), the degree-d certificates are excluded iff the quadric
system has only the trivial zero, i.e. the ideal is m-primary.  A modular m-primary test is rigorous: a nonzero
projective zero over Qbar reduces (after scaling to a primitive integral vector) to a nonzero zero mod any prime
above p for every p not dividing the clearing denominators, so m-primary mod p => only the trivial zero over C.
Requires rank_p J_{2d} = dim_Q J_{2d} (checked) so that the mod-p quadrics are reductions of the rational ones.

Usage: python3 b2_hessian.py d --cD yeven|sector:c1,c2|elem:g0,g1,g2,g3,g4,g5:v [--p P] [--kernel file.json ...]
       [--maxdeg 3] [--msolve] [--threads T] [--tag name]
"""
import sys, time, json, argparse, subprocess, os, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-excl')
from lib import *
from itertools import combinations_with_replacement
from math import comb

ap = argparse.ArgumentParser()
ap.add_argument('d', type=int)
ap.add_argument('--cD', default='yeven')
ap.add_argument('--p', type=int, default=P3)
ap.add_argument('--kernel', nargs='*', default=[])
ap.add_argument('--maxdeg', type=int, default=3)
ap.add_argument('--msolve', action='store_true')
ap.add_argument('--threads', type=int, default=2)
ap.add_argument('--tag', default='')
ap.add_argument('--nocheck3', action='store_true')
args = ap.parse_args()
d, p = args.d, args.p
t_start = time.time()
F, A, q, Fz, Az = load()
G = group_elements(Fz); gens = generators(G); chi = make_chi(gens); chiA = chi_of(Az, chi)
allchars = sorted(set(chi(m) for m in mons(6)))
assert len(allchars) == len(G) == 18
def g_of_char(c):
    # recover weight vector of an element from gens coefficients: not needed
    pass
# allowed D characters
def chi_g(m, g): return sum(a*b for a, b in zip(m, g)) % 6
yeven = [c for c in allchars if c[1] % 2 == 0]   # gens[1]^3 = tau
spec = args.cD
if spec == 'yeven':
    CD = set(yeven)
elif spec.startswith('sector:'):
    c = tuple(int(x) for x in spec[7:].split(',')); assert c in yeven; CD = {c}
elif spec.startswith('elem:'):
    body, v = spec[5:].rsplit(':', 1); g = tuple(int(x) for x in body.split(',')); v = int(v)
    # D is an eigenvector of g with chi_g(D) = v: find characters c with that property via a monomial representative
    rep = {}
    for m in mons(8):
        rep.setdefault(chi(m), m)
    CD = {c for c in yeven if chi_g(rep[c], g) == v}
else:
    raise SystemExit('bad --cD')
CD = sorted(CD)
print(f'd={d} p={p} |G|={len(G)} chiA={chiA} allowed D characters: {CD}', flush=True)

dF = [ideriv(Fz, i) for i in range(6)]
def chi_vf(i, m): return cadd(chi(m), cneg(chi(E6[i])))
# ---- unknowns ----
gauge = set()
for l in mons(d - 6):
    if chi(l) in CD:
        gauge.add((5, tuple(l[j] + (1 if j == 5 else 0) for j in range(6))))
wvars = [(i, m) for i in range(6) for m in mons(d - 5) if chi_vf(i, m) in CD and (i, m) not in gauge]
# kernel elements
kern = []
for fn in args.kernel:
    with open(fn) as f: K = json.load(f)
    items = K if isinstance(K, list) else [K]
    for it in items:
        Dk = core.parse(it['D']); Nk = [core.parse(s) for s in it['N']]
        # verify over Q: sum Nk_i F_i = A Dk
        lhs = {}
        for i in range(6): lhs = iadd(lhs, core.mul(Nk[i], core.pderiv(F, i)))
        assert core.sub(lhs, core.mul(A, Dk)) == {}, 'kernel element lift fails'
        # scale to integers: Dz = LD*Dk, Nz = Nk * LD*LA/LF, then common clearing
        LF = Fraction(next(iter(Fz.values()))) / Fraction(F[next(iter(Fz))])
        LA = Fraction(next(iter(Az.values()))) / Fraction(A[next(iter(Az))])
        LD = 1
        for c in Dk.values(): LD = lcm(LD, Fraction(c).denominator)
        Dz = {m: Fraction(c) * LD for m, c in Dk.items()}
        Nz = [{m: Fraction(c) * LD * LA / LF for m, c in n.items()} for n in Nk]
        L2 = 1
        for n in Nz:
            for c in n.values(): L2 = lcm(L2, c.denominator)
        for c in Dz.values(): L2 = lcm(L2, c.denominator)
        Dz = {m: int(c * L2) for m, c in Dz.items()}; Nz = [{m: int(c * L2) for m, c in n.items()} for n in Nz]
        lhs = {}
        for i in range(6): lhs = iadd(lhs, imul(Nz[i], dF[i]))
        assert iadd(lhs, imul(Az, Dz), -1) == {}, 'integral lift fails'
        cDk = chi_of(Dz, chi)
        if cDk in CD:
            kern.append((it.get('name', fn), Dz, Nz, cDk))
            print(f'  kernel element {it.get("name", fn)} char {cDk} terms {len(Dz)} lift terms {sum(len(n) for n in Nz)}', flush=True)
nK = len(kern); nW = len(wvars); n = nK + nW
weights = [k[3] for k in kern] + [chi_vf(i, m) for (i, m) in wvars]
print(f'unknowns: {nK} kernel + {nW} vector-field coordinates (gauge removed {len(gauge)}) = {n}', flush=True)

# ---- row sectors of degree 2d and left kernels of J_{2d} mod p ----
t0 = time.time()
rows_by = {}
for m in mons(2*d): rows_by.setdefault(chi(m), []).append(m)
Lmats = {}; rankJ_total = 0
for psi, rm in rows_by.items():
    ridx = {m: k for k, m in enumerate(rm)}
    cols = []
    for i in range(6):
        ci = cadd(psi, chi(E6[i]))
        for m in mons(2*d - 5):
            if chi(m) != ci: continue
            col = {}
            for e, v in dF[i].items():
                k = ridx[tuple(a+b for a, b in zip(m, e))]; col[k] = (col.get(k, 0) + v) % p
            cols.append(col)
    R, C = len(rm), len(cols)
    if C == 0:
        Lmats[psi] = (ridx, None, R); continue
    flat = [0] * (C * R)   # transpose: C x R, rows = J generators
    for j, col in enumerate(cols):
        for r, v in col.items(): flat[j*R + r] = v
    JT = nmod_mat(C, R, flat, p)
    X, nul = JT.nullspace()          # columns of X (R x nul) span {y : JT y = 0} = left kernel of J
    rankJ_total += R - nul
    Lmats[psi] = (ridx, X, nul)
ciJ = ci_dim_J(2*d)
print(f'J_{2*d}: rank mod p = {rankJ_total}, CI value = {ciJ}  ({round(time.time()-t0,1)}s)  dim R_{2*d} = {comb(2*d+5,5) - ciJ}', flush=True)
assert rankJ_total == ciJ, 'bad prime: Jacobian rank drop'

def reduce_vec(P):
    """Return (psi, list of residue coordinates mod p) of polynomial P (single sector) or None if P = 0."""
    if not P: return None
    psi = chi_of(P, chi); ridx, X, nul = Lmats[psi]
    if X is None: return (psi, [sum(0 for _ in [])])
    vec = [0] * len(ridx)
    for m, c in P.items(): vec[ridx[m]] = c % p
    out = []
    for k in range(nul):
        s = 0
        for r, v in enumerate(vec):
            if v: s += v * int(X[r, k])
        out.append(s % p)
    return (psi, out)

# ---- the quadratic form B_D(N) in (lam, w) ----
# variable a -> (D-part, N-part) as polys: for kernel k: (Dk, Nk); for w (i,m): (m*F_i, A*m*e_i)
Dparts = []; Nparts = []
for (_, Dz, Nz, _) in kern: Dparts.append(Dz); Nparts.append(Nz)
for (i, m) in wvars:
    Dparts.append(imul({m: 1}, dF[i]))
    Nv = [{} for _ in range(6)]; Nv[i] = imul(Az, {m: 1}); Nparts.append(Nv)
divN = []; 
for Nv in Nparts:
    s = {}
    for i in range(6): s = iadd(s, ideriv(Nv[i], i))
    divN.append(s)
dD = [[ideriv(Dp, j) for j in range(6)] for Dp in Dparts]
def Bcoef(a, b):
    """coefficient polynomial of x_a x_b (a <= b) in B_D(N) = D div N - N(D)."""
    P = imul(Dparts[a], divN[b])
    for j in range(6):
        if Nparts[b][j] and dD[a][j]: P = iadd(P, imul(Nparts[b][j], dD[a][j]), -1)
    if a != b:
        P = iadd(P, imul(Dparts[b], divN[a]))
        for j in range(6):
            if Nparts[a][j] and dD[b][j]: P = iadd(P, imul(Nparts[a][j], dD[b][j]), -1)
    return P
t0 = time.time()
quad = {}   # (psi, row) -> {(a,b): coeff}
npairs = 0
for a in range(n):
    for b in range(a, n):
        P = Bcoef(a, b)
        if not P: continue
        npairs += 1
        rv = reduce_vec(P)
        if rv is None: continue
        psi, out = rv
        for k, v in enumerate(out):
            if v: quad.setdefault((psi, k), {})[(a, b)] = v
quads = list(quad.values())
print(f'quadratic form: {npairs} nonzero coefficient polys, {len(quads)} nonzero quadrics mod p ({round(time.time()-t0,1)}s)', flush=True)
# weights of quadrics
def qweight(qd):
    (a, b) = next(iter(qd)); return cadd(weights[a], weights[b])
# sanity: gauge invariance -- random Euler gauge vector gives D = 6 l F; Phi must vanish on it (if gauge coords allowed)
# (we removed the gauge coordinates, so instead check that for random w the formula B_D(Aw) = -A w^T H w mod J holds)
rnd = random.Random(1)
wr = [rnd.randrange(p) for _ in range(n)]
# evaluate quadrics at wr: Phi(wr)
phi = [sum(c * wr[a] * wr[b] for (a, b), c in qd.items()) % p for qd in quads]
# independent evaluation: B_D(N) with D = sum x_a Dparts, N = sum x_a Nparts, computed directly mod p and reduced
Dtot = {}; Ntot = [{} for _ in range(6)]
for a in range(n):
    Dtot = iadd(Dtot, {m: c * wr[a] % p for m, c in Dparts[a].items()})
    for j in range(6): Ntot[j] = iadd(Ntot[j], {m: c * wr[a] % p for m, c in Nparts[a][j].items()})
Dtot = {m: c % p for m, c in Dtot.items() if c % p}; Ntot = [{m: c % p for m, c in Nj.items() if c % p} for Nj in Ntot]
divNt = {}
for j in range(6): divNt = iadd(divNt, ideriv(Ntot[j], j))
Bt = imul(Dtot, divNt, p)
for j in range(6): Bt = iadd(Bt, imul(Ntot[j], ideriv(Dtot, j), p), -1)
Bt = {m: c % p for m, c in Bt.items() if c % p}
# reduce Bt sector by sector
Bts = {}
for m, c in Bt.items(): Bts.setdefault(chi(m), {})[m] = c
ok_sanity = True
phi_direct = {}
for psi, P in Bts.items():
    rv = reduce_vec(P)
    for k, v in enumerate(rv[1]): phi_direct[(psi, k)] = v
for key, qd in quad.items():
    val = sum(c * wr[a] * wr[b] for (a, b), c in qd.items()) % p
    if val != phi_direct.get(key, 0): ok_sanity = False
for key, v in phi_direct.items():
    if v and key not in quad: ok_sanity = False
print(f'sanity (direct B_D(N) reduction at random point agrees with quadrics): {ok_sanity}', flush=True)
assert ok_sanity
# E1 sanity at the random point: sum N_i F_i = A D exactly mod p
lhs = {}
for j in range(6): lhs = iadd(lhs, imul(Ntot[j], dF[j], p))
e1ok = {m: c % p for m, c in iadd(lhs, imul(Az, Dtot, p), -1).items() if c % p} == {}
print(f'E1 at random point: {e1ok}', flush=True); assert e1ok

result = {'d': d, 'p': p, 'cD': spec, 'allowed_chars': CD, 'n_kernel': nK, 'n_w': nW, 'n_gauge_removed': len(gauge), 'n_vars': n,
          'n_quadrics': len(quads), 'rankJ_2d_mod_p': rankJ_total, 'ci_dimJ_2d': ciJ, 'sanity_direct': ok_sanity, 'kernel_names': [k[0] for k in kern]}
tag = args.tag or f'd{d}_{spec.replace(":", "_").replace(",", "")}'

# ---- degree-3 Macaulay test by weight blocks ----
if not args.nocheck3 and n > 0 and len(quads) > 0:
    t0 = time.time()
    # cubic monomials by weight
    cub = {}
    for a in range(n):
        for b in range(a, n):
            for c in range(b, n):
                cub.setdefault(cadd(weights[a], weights[b], weights[c]), []).append((a, b, c))
    full3 = True; blocks = {}
    for psi, cm in cub.items():
        cidx = {m: k for k, m in enumerate(cm)}
        cols = []
        for qd in quads:
            wq = qweight(qd)
            for v in range(n):
                if cadd(wq, weights[v]) != psi: continue
                col = {}
                for (a, b), c in qd.items():
                    key = tuple(sorted((a, b, v))); col[cidx[key]] = (col.get(cidx[key], 0) + c) % p
                cols.append(col)
        R, C = len(cm), len(cols)
        if C == 0: r = 0
        else:
            flat = [0] * (R * C)
            for j, col in enumerate(cols):
                for rr, v in col.items(): flat[rr*C + j] = v
            r = nmod_mat(R, C, flat, p).rank()
        blocks[str(psi)] = {'cubics': R, 'products': C, 'rank': r}
        if r < R: full3 = False
    result['degree3_blocks'] = blocks; result['m3_in_I'] = full3; result['t_deg3'] = round(time.time() - t0, 1)
    print(f'degree-3 Macaulay test: m^3 subset I = {full3}  blocks={blocks}  ({result["t_deg3"]}s)', flush=True)
# ---- msolve ----
if args.msolve and n > 0 and len(quads) > 0:
    fn = f'{PACKET}/msolve_{tag}.ms'
    names = [f'x{a}' for a in range(n)]
    with open(fn, 'w') as f:
        f.write(','.join(names) + '\n'); f.write(f'{p}\n')
        lines = []
        for qd in quads:
            terms = []
            for (a, b), c in qd.items():
                terms.append(f'{c}*{names[a]}*{names[b]}' if a != b else f'{c}*{names[a]}^2')
            lines.append('+'.join(terms))
        f.write(',\n'.join(lines) + '\n')
    t0 = time.time()
    outfn = f'{PACKET}/msolve_{tag}.out'
    cmd = ['msolve', '-v', '1', '-t', str(args.threads), '-g', '1', '-f', fn, '-o', outfn]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    # leading ideal: check pure powers for every variable
    lead = open(outfn).read()
    import re
    pure = set()
    for tok in re.findall(r'x(\d+)\^(\d+)|(?<![\w^])x(\d+)(?![\w^*])', lead.replace(' ', '')):
        if tok[0]: pure.add(int(tok[0]))
        elif tok[2]: pure.add(int(tok[2]))
    # tokens of form "xA^k" standing alone (no '*') are pure powers; need to parse lines properly
    pure = set()
    for mono in re.split(r'[,\n]', lead):
        mono = mono.strip().strip('[]')
        if not mono or '*' in mono: continue
        mm = re.fullmatch(r'x(\d+)(\^\d+)?', mono)
        if mm: pure.add(int(mm.group(1)))
    mprim = (len(pure) == n)
    result['msolve'] = {'seconds': round(time.time() - t0, 1), 'pure_powers': len(pure), 'n': n, 'm_primary': mprim, 'stderr_tail': pr.stderr[-2000:], 'stdout_tail': pr.stdout[-2000:]}
    print(f'msolve leading ideal: pure powers for {len(pure)}/{n} variables -> m-primary = {mprim} ({result["msolve"]["seconds"]}s)', flush=True)
result['seconds_total'] = round(time.time() - t_start, 1)
excluded = result.get('m3_in_I') or (result.get('msolve', {}).get('m_primary'))
result['excluded'] = bool(excluded)
receipt(f'hessian_{tag}', result)
log({'exp': 'hessian_relaxation', 'tag': tag, 'd': d, 'p': p, 'cD': spec, 'n_vars': n, 'n_quadrics': len(quads), 'm3_in_I': result.get('m3_in_I'),
     'msolve_m_primary': result.get('msolve', {}).get('m_primary'), 'excluded': result['excluded'], 'seconds': result['seconds_total']})
print('EXCLUDED' if excluded else 'NOT EXCLUDED (inconclusive)', result['seconds_total'], 's', flush=True)
