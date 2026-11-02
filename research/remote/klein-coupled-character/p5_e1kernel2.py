import sys, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc import *
F, A, q = load()
p = P0
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
a = 3
def e1_kernel(d, k, parD):
    Dm = wmonomials(d, parD, k)
    if not Dm: return 0, 0, 0
    pars = [(1 - parD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
    cols = [shift(Ap, m, 1, p) for m in Dm]
    for i in range(NV):
        kN = (a + k + W[i]) % 7
        for m in wmonomials(d + 1, pars[i], kN):
            cols.append(shift(Fi[i], m, p-1, p))
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0]*(R*C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e]*C + j] = v
    X, nullity = nmod_mat(R, C, flat, p).nullspace()
    nD = len(Dm)
    kerD = 0
    if nullity:
        kerD = nmod_mat(nD, nullity, [int(X[j, kk]) for j in range(nD) for kk in range(nullity)], p).rank()
    jd = 0
    if d >= 5:
        jc = []
        for i in range(NV):
            pm = (parD + (0 if i < 3 else 1)) % 2
            for m in wmonomials(d - 5, pm, (k + W[i]) % 7):
                jc.append(shift(Fi[i], m, 1, p))
        idx = {m: i for i, m in enumerate(Dm)}
        if jc:
            flat2 = [0]*(nD*len(jc))
            for j, col in enumerate(jc):
                for e, v in col.items():
                    if e in idx: flat2[idx[e]*len(jc)+j] = v
            jd = nmod_mat(nD, len(jc), flat2, p).rank()
    return nD, kerD, jd
for parD in (0, 1):
    print('--- D parity', 'Y-even' if parD == 0 else 'Y-odd')
    for d in range(1, 10):
        t = time.time(); line = []
        for k in range(7):
            nD, kerD, jd = e1_kernel(d, k, parD)
            line.append(f'k{k}:{kerD-jd}(+J{jd})/{nD}')
            log({'D': f'E1-kernel d={d} k={k} parD={parD}', 'd': d, 'family': 'e1-kernel', 'engine': 2, 'e1': None,
                 'consistent': None, 'rank': None, 'nullity': kerD - jd, 'seconds': round(time.time()-t, 2),
                 'note': f'nontrivial E1-kernel dim (mod J_d) = {kerD-jd}; J_d block dim {jd}; #monomials {nD}'})
        print(f'd={d}: ' + ' '.join(line) + f'  [{time.time()-t:.1f}s]', flush=True)
