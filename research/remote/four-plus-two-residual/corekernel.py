"""Core E1-kernel analysis: F = H(X0,X1,Y0,Y1) + B(X2,Y2); R(F) = R(H) (x) R(B).
For each core degree k (0..8) and Y-parity, compute dim ker( a : R(H)_k -> R(H)_{k+6} ), a = X0^3 Y0^3,
(as a subspace of the monomial span, modulo J(H)_k), with rational reconstruction via two primes."""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import monomials, shift, pderiv, reduce_mod, ratrecon, crt_pair
from flint import nmod_mat
from fractions import Fraction
F, A, q, rec = core.load_instance('four-plus-two-residual')
H = {e: c for e, c in F.items() if e[2] == 0 and e[5] == 0}
Bb = {e: c for e, c in F.items() if e[2] or e[5]}
assert core.add(H, Bb) == F
CORE = [0, 1, 3, 4]
def cmon(deg, par):
    return [e for e in monomials(deg, par) if e[2] == 0 and e[5] == 0]
def build(cols_polys, p):
    rows = {}
    for col in cols_polys:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols_polys)
    flat = [0] * (R * C)
    for j, col in enumerate(cols_polys):
        for e, v in col.items(): flat[rows[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p)
def kernel_basis(k, par, p):
    Hp = reduce_mod(H, p); ap = reduce_mod(A, p)
    Hi = [pderiv(Hp, i, p) for i in CORE]
    mons = cmon(k, par)
    tpar = (1 + par) % 2
    Jcols = []
    for idx, i in enumerate(CORE):
        for m2 in cmon(k + 1, (tpar - (0 if i < 3 else 1)) % 2):
            Jcols.append(shift(Hi[idx], m2, 1, p))
    Acols = [shift(ap, m, 1, p) for m in mons]
    MJ = build(Jcols, p); rJ = MJ.rank()
    MAJ = build(Acols + Jcols, p); rAJ = MAJ.rank()
    dimKer = len(mons) - (rAJ - rJ)
    # J(H)_k of parity par
    Jk = []
    if k >= 5:
        for idx, i in enumerate(CORE):
            for m3 in cmon(k - 5, (par - (0 if i < 3 else 1)) % 2):
                Jk.append(shift(Hi[idx], m3, 1, p))
    dimJk = build(Jk, p).rank() if Jk else 0
    # explicit kernel: nullspace of MAJ projected to A-coords, then rref
    X, nullity = MAJ.nullspace()
    nA = len(mons)
    vecs = [[int(X[j, kk]) for j in range(nA)] for kk in range(nullity)]
    vecs = [v for v in vecs if any(v)]
    basis = []
    if vecs:
        Mv = nmod_mat(len(vecs), nA, [x for v in vecs for x in v], p)
        Mr, rk = Mv.rref()
        for i in range(rk):
            basis.append({mons[j]: int(Mr[i, j]) for j in range(nA) if int(Mr[i, j])})
    return dict(k=k, par=par, nmons=len(mons), dimKer=dimKer, dimJk=dimJk, dimKer_modJ=dimKer - dimJk), basis, Jk
out = []
p1, p2 = core.PRIMES[0], core.PRIMES[1]
allbases = {}
for k in range(0, 9):
    for par in (0, 1):
        t0 = time.time()
        info1, b1, Jk1 = kernel_basis(k, par, p1)
        info2, b2, Jk2 = kernel_basis(k, par, p2)
        assert info1['dimKer'] == info2['dimKer']
        # rational reconstruction of rref basis (same pivots expected)
        rat = []
        ok = True
        for P1, P2 in zip(b1, b2):
            if set(P1) != set(P2): ok = False; break
            poly = {}
            for e in P1:
                a, m = crt_pair(P1[e], p1, P2[e], p2)
                r = ratrecon(a, m)
                if r is None: ok = False; break
                poly[e] = r
            if not ok: break
            rat.append(poly)
        info1['rational_lift'] = ok
        info1['seconds'] = round(time.time() - t0, 1)
        print(json.dumps(info1), flush=True)
        out.append(info1)
        if ok and rat:
            allbases[f'{k}_{par}'] = [core.to_str(P) for P in rat]
            if len(rat) <= 6:
                for P in rat: print('    ', core.to_str(P)[:200])
json.dump({'info': out, 'bases': allbases}, open('corekernel.json', 'w'), indent=1)
