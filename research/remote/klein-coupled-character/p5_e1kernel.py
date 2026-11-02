import sys, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc import *
F, A, q = load()
p = P0
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
a = 3
def e1_kernel(d, k, parD=0):
    """dim of {D deg d, parity parD, weight k : A D in J(F)} and dim J(F)_d in that block."""
    Dm = wmonomials(d, parD, k)
    pars = [(1 - parD - (0 if i < 3 else 1)) % 2 for i in range(NV)]  # A odd
    cols = []  # columns: D monomials (coefficient of A*m), then N_i monomials (-F_i * m)
    for m in Dm:
        cols.append(('D', shift(Ap, m, 1, p)))
    nN = 0
    for i in range(NV):
        kN = (a + k + W[i]) % 7
        for m in wmonomials(d + 1, pars[i], kN):
            cols.append(('N', shift(Fi[i], m, p-1, p))); nN += 1
    rows = {}
    for _, col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j, (_, col) in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    nD = len(Dm)
    if nullity == 0: kerD = 0
    else:
        proj = nmod_mat(nD, nullity, [int(X[j, kk]) for j in range(nD) for kk in range(nullity)], p)
        kerD = proj.rank()
    # dim J(F)_d in this block (trivial solutions)
    jd = 0
    if d >= 5:
        jc = []
        for i in range(NV):
            for m in wmonomials(d - 5, parD, (k + W[i]) % 7):
                jc.append(shift(Fi[i], m, 1, p))
        idx = {m: i for i, m in enumerate(Dm)}
        if jc:
            flat2 = [0]*(nD*len(jc))
            for j, col in enumerate(jc):
                for e, v in col.items():
                    if e in idx: flat2[idx[e]*len(jc)+j] = v
            jd = nmod_mat(nD, len(jc), flat2, p).rank()
    return nD, kerD, jd, nN, R
for d in range(1, 9):
    t = time.time()
    line = []
    for k in range(7):
        nD, kerD, jd, nN, R = e1_kernel(d, k)
        line.append(f'k{k}:{kerD}-{jd}/{nD}')
        log({'D': f'E1-kernel sweep d={d} k={k}', 'd': d, 'family': 'e1-kernel', 'e1': None, 'consistent': None,
             'rank': None, 'nullity': kerD - jd, 'seconds': round(time.time()-t, 2),
             'note': f'dim ker(A: R_d^(k) -> R_(d+6)) = {kerD} (incl. J_d dim {jd}) of {nD} even monomials'})
    print(f'd={d}: ' + '  '.join(line) + f'   [{time.time()-t:.1f}s]  (format kerdim-J_d dim / #monomials)', flush=True)
