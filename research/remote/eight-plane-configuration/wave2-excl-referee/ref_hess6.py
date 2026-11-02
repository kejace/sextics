"""Degree-6 Hessian relaxation rebuilt from scratch (own code) at a NEW prime, pure sectors and mixed sector.
w = linear vector field (w_i in S_1) with D = w(F) of parity par: w_{X_i} parity par, w_{Y_i} parity par+1.
Condition (necessary for a certificate, given E1-kernel_6 = J_6 and Ann_R(A)_6 = 0): Hess_F(w,w) in J_6 (even part in J_6^+, odd part in J_6^-).
Gauge: w -> w + c*Euler (D -> D + 6cF): slice = drop the Y2 d/dY2 coordinate. Quadrics = coordinates of NF_J(Hess(w,w)) over the standard
monomials of R_6^+ (and R_6^- for the mixed sector). msolve -g 1 -> leading ideal; m-primary (pure power of every variable) => V' = {0} over F_p
=> over Q by the Macaulay-rank argument (needs rank_p J_6 = 18 = exact, checked).  Singular dim/vdim as a second engine."""
import sys, time, re, subprocess, itertools, os
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl-referee')
import ref_lib as L
from flint import nmod_mat

p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000009
which = sys.argv[2:] or ['even', 'odd', 'mixed']
F, A, q = L.load()
Fp = L.red(F, p)
Fi = [L.pderiv(Fp, i, p) for i in range(L.NV)]
Fij = [[L.pderiv(Fi[i], j, p) for j in range(L.NV)] for i in range(L.NV)]
E1 = [tuple(1 if k == i else 0 for k in range(L.NV)) for i in range(L.NV)]

class Reducer:
    def __init__(self, par):
        self.mons = L.monomials(6, par); self.idx = {m: i for i, m in enumerate(self.mons)}
        gens = L.J_generators(Fi, 6, par, p)
        G = nmod_mat(len(gens), len(self.mons), p)
        for r, P in enumerate(gens):
            for e, v in P.items(): G[r, self.idx[e]] = v
        Gr, rk = G.rref(); self.rank = rk; self.exact = L.J_dim(6, par)
        piv = []
        for i in range(rk):
            j = 0
            while int(Gr[i, j]) == 0: j += 1
            piv.append(j)
        self.piv = piv; self.pivset = set(piv)
        self.std = [j for j in range(len(self.mons)) if j not in self.pivset]
        self.rows = [[int(Gr[i, j]) for j in self.std] for i in range(rk)]
    def nf(self, P):
        v = [0] * len(self.std); pc = {}
        sidx = {j: k for k, j in enumerate(self.std)}
        for e, c in P.items():
            j = self.idx[e]
            if j in self.pivset: pc[j] = (pc.get(j, 0) + c) % p
            else: v[sidx[j]] = (v[sidx[j]] + c) % p
        for i, j in enumerate(self.piv):
            c = pc.get(j, 0)
            if c:
                row = self.rows[i]
                for k in range(len(v)):
                    if row[k]: v[k] = (v[k] - c * row[k]) % p
        return v

def basis(par):
    B = []
    for i in range(L.NV):
        pm = (par + (0 if i < 3 else 1)) % 2
        for m in L.monomials(1, pm): B.append((i, m))
    return B

def hess(a, b):
    (i, m), (j, m2) = a, b
    return L.shift(Fij[i][j], tuple(x + y for x, y in zip(m, m2)), 1, p)

def write_ms(path, n, quadrics):
    with open(path, 'w') as f:
        f.write(','.join(f'w{a}' for a in range(n)) + '\n' + f'{p}\n')
        polys = []
        for qd in quadrics:
            terms = [(f'{c}*w{a}^2' if a == b else f'{c}*w{a}*w{b}') for (a, b), c in qd.items() if c]
            if terms: polys.append('+'.join(terms))
        f.write(',\n'.join(polys) + '\n')
    return len(polys)

def leading_ideal(path, n, threads=2):
    out = path.replace('.ms', '.lm'); lg = path.replace('.ms', '.msolve.log')
    t0 = time.time()
    rc = subprocess.run(f'msolve -v 1 -g 1 -t {threads} -f {path} -o {out} > {lg} 2>&1', shell=True).returncode
    lms = []
    for line in open(out):
        if line.startswith('#'): continue
        mon = [0] * n
        for a, k in re.findall(r'w(\d+)\^(\d+)', line): mon[int(a)] = int(k)
        if any(mon): lms.append(tuple(mon))
    pure = {a for lm in lms for a in range(n) if lm[a] and sum(lm) == lm[a]}
    # Hilbert function of S/LT(I): count standard monomials degree by degree (order ideal), stop when zero
    hf = []; cur = [tuple([0] * n)]; total = 0
    while cur and len(hf) < 40:
        hf.append(len(cur)); total += len(cur)
        nxt = set()
        for m in cur:
            for a in range(n):
                m2 = list(m); m2[a] += 1; m2 = tuple(m2)
                if m2 in nxt: continue
                if not any(all(m2[t] >= lm[t] for t in range(n)) for lm in lms): nxt.add(m2)
        cur = sorted(nxt)
    return {'rc': rc, 'n_leading': len(lms), 'vars_with_pure_power': len(pure), 'm_primary': len(pure) == n, 'hilbert_S_mod_LT': hf,
            'colength': total if not cur else None, 'msolve_seconds': round(time.time() - t0, 1)}

def singular_dim(path, n, timeout=1200):
    sing = path.replace('.ms', '.sing')
    lines = open(path).read().split('\n')
    polys = ''.join(lines[2:]).replace('\n', '')
    with open(sing, 'w') as f:
        f.write(f'ring r = {p}, (' + ','.join(f'w{a}' for a in range(n)) + '), dp;\n')
        f.write('ideal I = ' + polys.rstrip(',') + ';\n')
        f.write('ideal G = std(I);\nprint("dim"); print(dim(G)); print("vdim"); print(vdim(G)); quit;\n')
    t0 = time.time()
    try:
        r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=timeout)
        o = r.stdout.split()
        return {'singular_dim': int(o[o.index('dim') + 1]), 'singular_vdim': int(o[o.index('vdim') + 1]), 'singular_seconds': round(time.time() - t0, 1)}
    except Exception as ex:
        return {'singular_error': f'{type(ex).__name__}: {str(ex)[:120]}', 'singular_seconds': round(time.time() - t0, 1)}

red = {0: Reducer(0), 1: Reducer(1)}
L.log({'check': 'hess6: J_6 ranks at new prime', 'p': p, 'rank_J6_even': red[0].rank, 'exact_even': red[0].exact, 'rank_J6_odd': red[1].rank, 'exact_odd': red[1].exact,
       'dim_R6': [len(red[0].std), len(red[1].std)]})

for sector in which:
    t0 = time.time()
    if sector == 'mixed':
        B = basis(0) + basis(1)
    else:
        B = basis(0 if sector == 'even' else 1)
    gauge = [k for k, (i, m) in enumerate(B) if i == 5 and m == E1[5]]      # Y2 d/dY2 coordinate (present in the even sector only)
    keep = [k for k in range(len(B)) if k not in gauge]
    n = len(keep)
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    quadrics_by_col = {}
    for (a, b) in pairs:
        h = hess(B[keep[a]], B[keep[b]])
        if a != b: h = {e: (2 * v) % p for e, v in h.items()}
        parts = {0: {}, 1: {}}
        for e, v in h.items(): parts[L.ydeg(e) % 2][e] = v
        for par in (0, 1):
            if not parts[par]: continue
            if sector != 'mixed' and par == 1: raise RuntimeError('odd part in a pure sector')
            v = red[par].nf(parts[par])
            for k, c in enumerate(v):
                if c: quadrics_by_col.setdefault((par, k), {})[(a, b)] = c
    quadrics = list(quadrics_by_col.values())
    # linearization rank = rank of the (pairs x quadrics) coefficient matrix
    Mq = nmod_mat(len(pairs), len(quadrics), p); pidx = {pr: i for i, pr in enumerate(pairs)}
    for j, qd in enumerate(quadrics):
        for pr, c in qd.items(): Mq[pidx[pr], j] = c
    rk = Mq.rank()
    path = f'{L.REF}/ref_hess6_{sector}_p{p}.ms'
    nq = write_ms(path, n, quadrics)
    info = {'check': f'hess6 {sector} sector, own build', 'p': p, 'dimW': len(B), 'gauge_removed': len(gauge), 'n_vars': n, 'dimSym2': len(pairs),
            'n_quadrics_written': nq, 'rank_linearization': rk, 'linearization_proves_empty': rk == len(pairs), 't_build': round(time.time() - t0, 1), 'msolve_input': path}
    info.update(leading_ideal(path, n))
    info.update(singular_dim(path, n))
    info['verdict'] = ("V' = {0}: no coprime degree-6 denominator in this sector (rigorous over C by the Macaulay-rank transfer, J_6 no-drop checked)"
                       if info['m_primary'] and red[0].rank == red[0].exact and red[1].rank == red[1].exact else 'NOT m-primary / check failed')
    L.log(info)
