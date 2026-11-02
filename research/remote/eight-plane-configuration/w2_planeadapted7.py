"""Hessian quadrics at d=7 restricted to the plane-adapted subspace: w = sum_j (L_j R_j + M_j S_j) + sum_j Q_j v_j."""
import sys, json, subprocess, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, mul, add, scale
from fractions import Fraction
WD = '/tmp/hodge/results/eight-plane-configuration'
P = core.PRIMES[1]
def lin(i):
    e = [0]*6; e[i] = 1; return {tuple(e): Fraction(1)}
Z = [lin(i) for i in range(6)]
Q = [add(mul(Z[i], Z[i]), mul(Z[i+3], Z[i+3])) for i in range(3)]
S = [[Z[j] if j in (i, i+3) else {} for j in range(6)] for i in range(3)]
R = [[(scale(Z[i+3], -1) if j == i else (Z[i] if j == i+3 else {})) for j in range(6)] for i in range(3)]
for par in (0, 1):
    dat = json.load(open(f'{WD}/hess7local-quadrics-d7-p{par}.json'))
    Wm = [(i, tuple(m)) for i, m in dat['Wm']]; widx = {wm: k for k, wm in enumerate(Wm)}
    rows = [{(a, b): v for (a, b), v in r} for r in dat['rows']]
    # generators of the subspace: list of (name, field comps)
    gens = []
    Ls = [Z[3], Z[4], Z[5]] if par == 0 else [Z[0], Z[1], Z[2]]   # L_j parity: even D -> L odd
    Ms = [Z[0], Z[1], Z[2]] if par == 0 else [Z[3], Z[4], Z[5]]
    for j in range(3):
        for L in Ls: gens.append((f'L*R{j}', [mul(L, c) if c else {} for c in R[j]]))
        for M in Ms: gens.append((f'M*S{j}', [mul(M, c) if c else {} for c in S[j]]))
        for b in (range(3) if par == 0 else range(3, 6)): gens.append((f'Q{j}*d{b}', [Q[j] if i == b else {} for i in range(6)]))
    # coordinate vectors mod p
    T = []
    for name, comps in gens:
        vec = {}
        for i in range(6):
            for m, c in comps[i].items():
                cc = Fraction(c); vec[widx[(i, m)]] = (vec.get(widx[(i, m)], 0) + cc.numerator * pow(cc.denominator, -1, P)) % P
        T.append(vec)
    n = len(T)
    # substitute: quadric q(w) = sum Q[kl] w_k w_l with w = sum_s t_s T_s  -> coefficients in t
    polys = []
    for row in rows:
        coef = {}
        for (k, l), v in row.items():
            for s in range(n):
                ak = T[s].get(k, 0)
                if not ak: continue
                for t in range(n):
                    al = T[t].get(l, 0)
                    if not al: continue
                    key = (min(s, t), max(s, t))
                    coef[key] = (coef.get(key, 0) + v * ak * al) % P
        coef = {k: v for k, v in coef.items() if v}
        if coef: polys.append(coef)
    # rank of the restricted quadric system
    keys = sorted({k for c in polys for k in c}); kidx = {k: i for i, k in enumerate(keys)}
    from flint import nmod_mat
    flat = [0] * (len(polys) * len(keys))
    for r, c in enumerate(polys):
        for k, v in c.items(): flat[r * len(keys) + kidx[k]] = v
    rk = nmod_mat(len(polys), len(keys), flat, P).rank() if polys else 0
    names = ','.join(f't{s}' for s in range(n))
    terms = []
    for c in polys:
        terms.append('+'.join(f'{v}*t{s}*t{t}' if s != t else f'{v}*t{s}^2' for (s, t), v in c.items()))
    sing = f'''LIB "primdec.lib";
ring r = {P}, ({names}), dp;
option(redSB);
ideal I = {','.join(terms)};
ideal G = std(I);
"DIM", dim(G), "DEGREE", degree(G);
int k; for (k = 0; k <= 6; k++) {{ "HF", k, size(kbase(G, k)); }}
list L = primdecGTZ(G);
"NCOMP", size(L);
int i; for (i = 1; i <= size(L); i++) {{ ideal Pi = std(L[i][2]); "COMP", i, "dim", dim(Pi), "deg", degree(Pi); Pi; }}
quit;
'''
    open(f'{WD}/planeadapted7-p{par}.sing', 'w').write(sing)
    out = subprocess.run(['Singular', '-q', f'{WD}/planeadapted7-p{par}.sing'], capture_output=True, text=True, timeout=600).stdout
    print(f'par {par}: subspace dim {n}, restricted quadrics {len(polys)} (rank {rk})')
    print('\n'.join(l for l in out.splitlines() if l.startswith(('DIM', 'HF', 'NCOMP', 'COMP', 'Pi[', '_['))))
    with open(f'{WD}/wave2-search.log.jsonl', 'a') as f:
        f.write(json.dumps({'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'd=7 Hessian quadrics restricted to plane-adapted fields (L*R_j, M*S_j, Q_j*const)', 'd': 7, 'par': par,
                            'subspace_dim': n, 'n_quadrics': len(polys), 'rank': rk, 'singular': [l for l in out.splitlines() if l.startswith(('DIM', 'NCOMP', 'COMP'))]}) + '\n')
