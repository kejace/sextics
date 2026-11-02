"""Construction-lens helpers (wave 2): restricted E1 kernels on structured denominator families."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
from hodge.core import monomials, shift, reduce_mod
from flint import nmod_mat
from fractions import Fraction

IID = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
c = Fraction(5, 7**7)
L = core.parse('X0 + 2*X1 + 3*X2'); M = core.parse('Y0 + 3*Y1 + 5*Y2')
eX = core.parse('X0**6 + X1**6 + X2**6'); eY = core.parse('Y0**6 + Y1**6 + Y2**6')

def ppow(P, k):
    return core.power(P, k) if k > 0 else {(0,)*6: Fraction(1)}

def fam_vert(d, par):
    """monomials e(X)^a e(Y)^b L^i M^j, 6a+6b+i+j=d, j%2==par  (pullbacks of curves on the (eps_X,eps_Y,t) base)."""
    out = []
    for a in range(d // 6 + 1):
        for b in range(d // 6 + 1):
            r = d - 6 * a - 6 * b
            if r < 0: continue
            for j in range(r + 1):
                i = r - j
                if j % 2 != par: continue
                P = core.mul(core.mul(ppow(eX, a), ppow(eY, b)), core.mul(ppow(L, i), ppow(M, j)))
                out.append((f'eX^{a} eY^{b} L^{i} M^{j}', P))
    return out

def fam_coneXM(d, par):
    """C[X0,X1,X2,M]_d with M-degree parity par (cones with vertex the line X=0,M=0; pullbacks of surfaces from the (x,t) base)."""
    out = []
    for j in range(par, d + 1, 2):
        Mj = ppow(M, j)
        for m in monomials(d - j, 0):
            if any(m[3:]): continue
            out.append((f'X^{m[:3]} M^{j}', shift(Mj, m, 1)))
    return out

def fam_coneLY(d, par):
    """C[L,Y0,Y1,Y2]_d with Y-parity par (cones with vertex the line Y=0,L=0; pullbacks of surfaces from the (y,t) base)."""
    out = []
    for i in range(0, d + 1):
        Li = ppow(L, i)
        for m in monomials(d - i, par):
            if any(m[:3]): continue
            out.append((f'L^{i} Y^{m[3:]}', shift(Li, m, 1)))
    return out

def restricted_kernel(F, A, fam, d, p, verbose=True):
    """{D in span(fam) : A*D in J(F)_{d+6}} mod p. Returns (kernel vectors as coefficient dicts over fam names, info)."""
    t0 = time.time()
    Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
    parA = core.parity_of(Ap)
    names = [n for n, _ in fam]
    polys = [reduce_mod(P, p) for _, P in fam]
    pars = {core.parity_of(P) for P in polys if P}
    assert len(pars) == 1, pars
    parD = pars.pop()
    target_par = (parA + parD) % 2
    gens = jacobian.jacobian_generators(Fp, d + 6, p, target_par)
    cols = [core.mul(Ap, P, p) for P in polys] + [g[2] for g in gens]
    rows = {}
    for col in cols:
        for e in col: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for e, v in col.items(): flat[rows[e] * C + j] = v
    Mx = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nul = Mx.nullspace()
    t2 = time.time()
    K = len(fam)
    vecs = []
    if nul:
        proj = nmod_mat(nul, K, [int(X[j, k]) for k in range(nul) for j in range(K)], p)
        Pr, rk = proj.rref()
        for i in range(rk):
            vecs.append({names[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])})
    info = {'d': d, 'parD': parD, 'fam_dim': K, 'rows': R, 'cols': C, 'nullity_total': nul,
            'kernel_dim': len(vecs), 'p': p, 't_build': round(t1 - t0, 1), 't_null': round(t2 - t1, 1)}
    if verbose: print(json.dumps(info), flush=True)
    return vecs, info

def vec_to_poly(vec, fam, p):
    name2P = {n: reduce_mod(P, p) for n, P in fam}
    out = {}
    for n, cf in vec.items():
        out = core.add(out, core.scale(name2P[n], cf, p), p)
    return out

def analyze_kernel(F, vecs, fam, d, p):
    """For each kernel vector: is D a multiple of F? is D in J(F)_d? (both mod p)."""
    res = []
    Fp = reduce_mod(F, p)
    for v in vecs:
        D = vec_to_poly(v, fam, p)
        divF, _ = jacobian.divisible_by(Fp, D, p) if d >= 6 else (False, None)
        inJ, _ = jacobian.in_jacobian(Fp, D, p) if d >= 5 else (False, None)
        res.append({'vec': {k: (x if x <= p // 2 else x - p) for k, x in v.items()}, 'F_divides': bool(divF), 'in_J': bool(inJ), 'nterms': len(D)})
    return res
