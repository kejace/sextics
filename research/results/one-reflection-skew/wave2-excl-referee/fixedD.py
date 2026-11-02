"""Independent fixed-D exact test (6-var, sector-restricted) and strict 5-var test; exact ranks over Z."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl-referee')
from reflib import *
from flint import nmod_mat, fmpz_mat
LOG = '/tmp/hodge/results/one-reflection-skew/wave2-excl-referee/referee.log.jsonl'
UNIT = [tuple(1 if k==i else 0 for k in range(NV)) for i in range(NV)]

def system6(D):
    d = pdeg(D); chD = pchi(D); assert chD is not None
    Ncols = [(i, m) for i in range(NV) for m in secmons(d+1, chi_add(CHI_A, chD, chi(UNIT[i])))]
    Qcols = secmons(2*d-6, chi_add(CHI_A, chD, chD))
    r1 = secmons(d+6, chi_add(CHI_A, chD)); r2 = secmons(2*d, chi_add(CHI_A, chD, chD))
    i1 = {e:k for k,e in enumerate(r1)}; i2 = {e:len(r1)+k for k,e in enumerate(r2)}
    Di = [pderiv(D, i) for i in range(NV)]
    cols = []
    for (i, m) in Ncols:
        mm = mono(m); col = {}
        for e, v in pmul(mm, Fi[i]).items(): col[i1[e]] = v                       # E1: N_i F_i
        P2 = padd(pmul(D, pderiv(mm, i)), pmul(mm, Di[i]), -1)                     # E2: D d_i(m) - m D_i
        for e, v in P2.items(): col[i2[e]] = col.get(i2[e],0) + v
        cols.append(col)
    for m in Qcols:
        col = {}
        for e, v in pmul(F7, mono(m)).items(): col[i2[e]] = -v                   # E2: - F Q
        cols.append(col)
    nr = len(r1) + len(r2)
    b = {}
    for e, v in pmul(A, D).items(): b[i1[e]] = v
    return cols, b, nr, dict(nN=len(Ncols), nQ=len(Qcols), rowsE1=len(r1), rowsE2=len(r2), d=d, chi=chD)

def system5(D):
    d = pdeg(D); chD = pchi(D)
    Ncols = [(i, m) for i in range(5) for m in secmons(d+1, chi_add(CHI_A, chD, chi(UNIT[i])), 5)]
    r1 = secmons(d+6, chi_add(CHI_A, chD), 5); r2 = secmons(2*d, chi_add(CHI_A, chD, chD), 5)
    i1 = {e:k for k,e in enumerate(r1)}; i2 = {e:len(r1)+k for k,e in enumerate(r2)}
    Di = [pderiv(D, i) for i in range(5)]
    cols = []
    for (i, m) in Ncols:
        mm = mono(m); col = {}
        for e, v in pmul(mm, fi[i]).items(): col[i1[e]] = v
        P2 = padd(pmul(D, pderiv(mm, i)), pmul(mm, Di[i]), -1)
        for e, v in P2.items(): col[i2[e]] = col.get(i2[e],0) + v
        cols.append(col)
    b = {}
    for e, v in pmul(A, D).items(): b[i1[e]] = v
    return cols, b, len(r1)+len(r2), dict(nN=len(Ncols), nQ=0, rowsE1=len(r1), rowsE2=len(r2), d=d, chi=chD)

def to_mats(cols, b, nr, p=None):
    nc = len(cols)
    M = fmpz_mat(nr, nc) if p is None else nmod_mat(nr, nc, p)
    Mb = fmpz_mat(nr, nc+1) if p is None else nmod_mat(nr, nc+1, p)
    for j, col in enumerate(cols):
        for r, v in col.items():
            vv = v if p is None else v % p
            M[r, j] = vv; Mb[r, j] = vv
    for r, v in b.items(): Mb[r, nc] = v if p is None else v % p
    return M, Mb

def run(name, D, five=False, exact=True, primes=(1000000007, 998244353)):
    t0 = time.time()
    cols, b, nr, meta = (system5 if five else system6)(D)
    out = dict(exp='referee-fixed-D' + ('-5var' if five else ''), D=name, **meta, rows=nr, cols=len(cols))
    for p in primes:
        M, Mb = to_mats(cols, b, nr, p)
        out[f'rank_p{p}'] = (M.rank(), Mb.rank())
    print(name, {k:v for k,v in out.items() if k.startswith('rank') or k in ('rows','cols','nN','nQ')}, f'{time.time()-t0:.1f}s', flush=True)
    if exact:
        t1 = time.time()
        M, Mb = to_mats(cols, b, nr, None)
        rM = M.rank(); rMb = Mb.rank()
        out['rank_exact'] = (rM, rMb); out['exact_seconds'] = round(time.time()-t1, 1)
        out['verdict'] = 'EXCLUDED (exact over Z: rank[M|b] = rank M + 1)' if rMb == rM + 1 else ('consistent over Q' if rMb == rM else '??')
        print('  exact:', out['rank_exact'], out['verdict'], f'{time.time()-t1:.1f}s', flush=True)
    out['seconds'] = round(time.time()-t0, 1)
    log(LOG, **out)
    return out

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'Y1Y2^6'
    Y2_6 = mono((0,0,0,0,0,6))
    cands = {'Y1Y2^6': mono((0,0,0,0,1,6)), 'Y0Y2^6': mono((0,0,0,1,0,6)), 'X1Y2^6': mono((0,1,0,0,0,6)),
             'Y1f': pmul(mono((0,0,0,0,1,0)), f7), 'Y0f': pmul(mono((0,0,0,1,0,0)), f7), 'X1f': pmul(mono((0,1,0,0,0,0)), f7)}
    exact = '--noexact' not in sys.argv
    for w in which.split(','):
        run(w, cands[w], five=w.endswith('f'), exact=exact)
