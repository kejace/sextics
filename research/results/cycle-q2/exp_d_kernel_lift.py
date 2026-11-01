"""Experiment D: exact rational E1-kernel at degree d (per character sector), reduced modulo J_d,
then full fixed-D tests (D and D^2) for every complement representative and small pencils."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np, math
from hodge import core
from fractions import Fraction
from flint import fmpz_mat, fmpq_mat, fmpq

d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
DO_FULL = (len(sys.argv) <= 2) or sys.argv[2] != 'nofull'
F, A, q, rec = core.load_instance('cycle-q2')
Gw = eq.symmetry_group(F); gens = eq.generators(Gw); ng = gens.shape[0]
kA = eq.key_of_poly(A, gens)
kx = [int(eq.keys_of(np.eye(6, dtype=np.int64)[i:i+1], gens)[0]) for i in range(6)]
F2 = core.scale(F, Fraction(2))            # integer coefficients
Fi = [core.pderiv(F2, i) for i in range(6)]  # integer polys
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']

def exps_in_sector(deg, kappa):
    allm = eq.monomial_exps(deg); km = eq.keys_of(allm, gens)
    return [tuple(int(x) for x in e) for e in allm[km == kappa]]

def J_sector(deg, kappa):
    """list of integer polys x_j F_i (deg) with character kappa"""
    out = []
    for i in range(6):
        need = eq.key_add(kappa, kx[i], ng)           # key(m) = kappa + key(x_i)
        for m in exps_in_sector(deg - 5, need):
            out.append(core.shift(Fi[i], m, 1))
    return out

def poly_vec(P, idx):
    v = [0] * len(idx)
    for e, c in P.items():
        v[idx[e]] = int(c)
    return v

def rational_kernel(d, kappa):
    """exact: {D in (S_d)_kappa : A D in J}; returns (list of rational D polys spanning K), dim"""
    Dm = exps_in_sector(d, kappa)
    if not Dm: return [], 0
    cols = [core.mul(A, {m: Fraction(1)}) for m in Dm]
    meta = len(Dm)
    for i in range(6):
        need = eq.key_add(eq.key_add(kA, kappa, ng), kx[i], ng)
        for m in exps_in_sector(d + 7, need):
            cols.append(core.shift(Fi[i], m, 1))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    M = [[0] * len(cols) for _ in range(len(rows))]
    for j, c in enumerate(cols):
        for e, v in c.items(): M[rows[e]][j] = int(v)
    Z = fmpz_mat(M)
    X, nul = Z.nullspace()
    # D-components
    vecs = []
    for kcol in range(nul):
        v = [int(X[j, kcol]) for j in range(meta)]
        if any(v): vecs.append(v)
    if not vecs: return [], 0
    R = fmpq_mat([[fmpq(x) for x in v] for v in vecs])
    Rr, rk = R.rref()
    basis = []
    for i in range(rk):
        P = {}
        for j in range(meta):
            c = Rr[i, j]
            if c != 0: P[Dm[j]] = Fraction(int(c.p), int(c.q))
        basis.append(P)
    return basis, rk

def complement_mod_J(basis, d, kappa):
    """representatives of K / (K cap J) : rref of [J ; K] and keep K-rows that add rank"""
    Dm = exps_in_sector(d, kappa); idx = {m: j for j, m in enumerate(Dm)}
    J = J_sector(d, kappa)
    Jv = [poly_vec(P, idx) for P in J]
    rJ = fmpq_mat([[fmpq(x) for x in v] for v in Jv]).rank() if Jv else 0
    reps = []
    cur = list(Jv)
    currank = rJ
    for P in basis:
        v = [0] * len(Dm)
        for e, c in P.items(): v[idx[e]] = c
        # scale to integers
        den = 1
        for c in v:
            if c: den = math.lcm(den, Fraction(c).denominator)
        vi = [int(Fraction(c) * den) for c in v]
        test = cur + [vi]
        r = fmpq_mat([[fmpq(x) for x in w] for w in test]).rank()
        if r > currank:
            cur = test; currank = r
            reps.append({Dm[j]: Fraction(vi[j]) for j in range(len(Dm)) if vi[j]})
    return reps, rJ

def clean(P):
    # make coefficients integers with gcd 1
    den = 1
    for c in P.values():
        den = math.lcm(den, c.denominator)
    g = 0
    for c in P.values():
        g = math.gcd(g, int(c * den))
    return {e: Fraction(int(c * den) // g) for e, c in P.items()}

t0 = time.time()
allD = eq.monomial_exps(d); kD = eq.keys_of(allD, gens)
keys = sorted(set(kD.tolist()))
summary = []
cands = []
for kappa in keys:
    basis, rk = rational_kernel(d, kappa)
    if rk == 0: continue
    reps, rJ = complement_mod_J(basis, d, kappa)
    if not reps: continue
    reps = [clean(P) for P in reps]
    # LLL-reduce the complement representatives together with J to get small coefficients
    Dm = exps_in_sector(d, kappa); idx = {m: j for j, m in enumerate(Dm)}
    L = fmpz_mat([poly_vec(clean(P), idx) for P in reps]).lll()
    reps = [clean({Dm[j]: Fraction(int(L[i, j])) for j in range(len(Dm)) if int(L[i, j])}) for i in range(len(reps))]
    print(f'key {kappa}: dim K={rk} dim J={rJ} complement {len(reps)}: ' + ' | '.join(core.to_str(P)[:120] for P in reps), flush=True)
    summary.append({'key': kappa, 'dimK': rk, 'dimJ': rJ, 'reps': [core.to_str(P) for P in reps]})
    for P in reps: cands.append((kappa, P, 'rep'))
    if len(reps) == 2:
        for t in [1, -1, 2, -2]:
            cands.append((kappa, core.add(reps[0], core.scale(reps[1], Fraction(t))), f'pencil t={t}'))
    if len(reps) >= 3:
        rng = np.random.default_rng(kappa)
        for _ in range(3):
            c = rng.integers(-3, 4, size=len(reps))
            P = {}
            for ci, Pi in zip(c, reps):
                if ci: P = core.add(P, core.scale(Pi, Fraction(int(ci))))
            if P: cands.append((kappa, P, f'random combo {c.tolist()}'))
json.dump(summary, open(f'/tmp/hodge/results/cycle-q2/exp_d_kernel_d{d}.json', 'w'), indent=1)
print(f'total candidates: {len(cands)}  [{time.time()-t0:.0f}s]', flush=True)
if not DO_FULL: sys.exit(0)
for kappa, P, how in cands:
    Ds = core.to_str(P)
    r = eq.test_D_fast(F, A, q, P, f'E1-kernel-d{d}', note=f'key {kappa} {how} [D]', verbose=False, do_e1=True)
    s = f'  D: cons={r.get("consistent")} rk={r.get("rank")} nul={r.get("nullity")} cols={r.get("cols")} |H|={r.get("order")} ({r["seconds"]}s)'
    if r.get('consistent'):
        print('!!!!! CONSISTENT (D) !!!!!', Ds, flush=True)
    if r.get('e1') is not False:
        rb = eq.test_D_fast(F, A, q, P, f'E1-kernel-d{d}', note=f'key {kappa} {how} [D] 2nd prime', verbose=False, do_e1=False, p=4194287)
        s += f' (p2: cons={rb.get("consistent")} rk={rb.get("rank")})'
        if rb.get('consistent'): print('!!!!! CONSISTENT (D, p2) !!!!!', Ds, flush=True)
    if d <= 6:
        P2 = core.mul(P, P)
        r2 = eq.test_D_fast(F, A, q, P2, f'E1-kernel-d{d}', note=f'key {kappa} {how} [D^2]', verbose=False, Dstr='(' + Ds + ')**2')
        s += f' | D^2: cons={r2.get("consistent")} rk={r2.get("rank")} nul={r2.get("nullity")} cols={r2.get("cols")} ({r2["seconds"]}s)'
        if r2.get('consistent'):
            print('!!!!! CONSISTENT (D^2) !!!!!', Ds, flush=True)
    print(f'key {kappa} {how}: {Ds[:90]}' + s, flush=True)
print('DONE', flush=True)
