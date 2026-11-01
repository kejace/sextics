"""E1 kernel sweep: dim{D Y-even, deg d : A*D in J(F)} and dim J_d^even, mod p.
Also saves a basis (mod p) of the kernel modulo J_d as sparse dicts (for later inspection)."""
import sys, time, json, pickle
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
from hodge import core
from hodge.core import monomials, mul, shift, pderiv
from flint import nmod_mat
from orslib import F, A, q, log
p = core.PRIMES[0]
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]

def e1_kernel(d, verbose=True, save=True):
    t0 = time.time()
    Dmons = monomials(d, 0)
    Ncols = []
    for i in range(6):
        par = 1 if i < 3 else 0
        for m in monomials(d + 1, par):
            Ncols.append((i, m))
    rows = {}
    colpolys = []
    for m in Dmons:
        colpolys.append(shift(Ap, m, 1, p))
    for (i, m) in Ncols:
        colpolys.append(shift(Fi[i], m, 1, p))
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(colpolys)
    nD = len(Dmons)
    flat = [0] * (R * C)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    # projection to D coords
    if nullity == 0:
        kerD = 0; basisD = []
    else:
        proj = nmod_mat(nD, nullity, [int(X[j, k]) for j in range(nD) for k in range(nullity)], p)
        T = proj.transpose()
        Tr, kerD = T.rref()
        basisD = []
        for i in range(kerD):
            poly = {}
            for j in range(nD):
                v = int(Tr[i, j])
                if v: poly[Dmons[j]] = v
            basisD.append(poly)
    # dim J_d^even: span of m*F_i with m of degree d-5 and parity so that product is even
    Jcols = []
    if d >= 5:
        for i in range(6):
            par = 0 if i < 3 else 1
            for m in monomials(d - 5, par):
                Jcols.append(shift(Fi[i], m, 1, p))
    if Jcols:
        rowsJ = {e: j for j, e in enumerate(Dmons)}
        flatJ = [0] * (nD * len(Jcols))
        for j, col in enumerate(Jcols):
            for e, v in col.items():
                flatJ[rowsJ[e] * len(Jcols) + j] = v
        MJ = nmod_mat(nD, len(Jcols), flatJ, p)
        dimJ = MJ.rank()
    else:
        dimJ = 0
    # reduce kernel basis modulo J: rank of [J | ker] minus dimJ gives kernel mod J
    kmodJ = kerD - dimJ  # since J_d subset kernel
    t3 = time.time()
    out = {'d': d, 'rows': R, 'cols': C, 'nD': nD, 'nullity_full': nullity, 'kerD': kerD,
           'dimJ_even': dimJ, 'ker_mod_J': kmodJ, 't_build': round(t1 - t0, 1),
           't_null': round(t2 - t1, 1), 't_total': round(t3 - t0, 1)}
    if verbose:
        print(json.dumps(out), flush=True)
    if save:
        with open(f'/tmp/hodge/results/one-reflection-skew/e1ker_d{d}.pkl', 'wb') as f:
            pickle.dump({'p': p, 'basisD': basisD, 'Dmons': Dmons, 'info': out}, f)
    log(f'<E1 kernel sweep degree {d}>', d, 'e1-kernel-sweep', kerD > 0, None, kerD, kmodJ, t3 - t0,
        f'dim{{D even deg {d}: AD in J}}={kerD}, dim J_{d}^even={dimJ}, kernel mod J={kmodJ}, matrix {R}x{C}')
    return out

if __name__ == '__main__':
    ds = [int(a) for a in sys.argv[1:]] or list(range(1, 9))
    for d in ds:
        e1_kernel(d)
