"""Structured denominator families at a given degree, mod p (p = 1 mod 4 so that i exists)."""
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *
from hodge.core import monomials, mul, add, sub, scale, power, reduce_mod, parity_of
from flint import nmod_mat

def sqrt_m1(p):
    for a in range(2, p):
        r = pow(a, (p-1)//4, p)
        if r*r % p == p-1: return r
    raise ValueError

def Kq():
    d = json.load(open(f'{WD}/K_decomp.json'))
    return [core.parse(s) for s in d['K']], [core.parse(s) for s in d['H']]

def plane_forms(eps, p):
    """u_j = Y_j - eps_j i X_j, ubar_j = Y_j + eps_j i X_j mod p."""
    I = sqrt_m1(p)
    u, ub = [], []
    for j in range(3):
        Xj = {tuple(1 if k==j else 0 for k in range(6)): 1}
        Yj = {tuple(1 if k==j+3 else 0 for k in range(6)): 1}
        u.append(add(Yj, scale(Xj, (-eps[j]*I) % p, p), p))
        ub.append(add(Yj, scale(Xj, (eps[j]*I) % p, p), p))
    return u, ub

def poly_pow_monos(forms, deg, p):
    """all monomials of degree deg in the given forms (list of polys), mod p"""
    out = {}
    n = len(forms)
    for a in itertools.product(range(deg+1), repeat=n):
        if sum(a) != deg: continue
        P = {tuple([0]*6): 1}
        for f, k in zip(forms, a):
            if k: P = mul(P, power(f, k, p), p)
        out[a] = P
    return out

def families(d, par, p):
    """returns dict name -> list of polys (mod p) of degree d and Y-parity par"""
    K, H = Kq()
    Kp = [reduce_mod(k, p) for k in K]; Hp = [reduce_mod(h, p) for h in H]
    diag = [reduce_mod(P(f'X{i}**4 - X{i}**2*Y{i}**2 + Y{i}**4'), p) for i in range(3)]
    Qp = [reduce_mod(Q0, p), reduce_mod(Q1, p), reduce_mod(Q2, p)]
    fams = {}
    if par == 0 and d % 2 == 0:
        fams['G_k(Q)'] = list(poly_pow_monos(Qp, d//2, p).values())
        # subring C[Q,K,H,diag]: Q^alpha * (product of quartics)
        quartics = Kp + Hp + diag
        sub_ = []
        for nq in range(0, d//4 + 1):
            dq = d - 4*nq
            if dq < 0 or dq % 2: continue
            Qm = poly_pow_monos(Qp, dq//2, p)
            for combo in itertools.combinations_with_replacement(range(9), nq):
                Pq = {tuple([0]*6): 1}
                for c in combo: Pq = mul(Pq, quartics[c], p)
                for qm in Qm.values(): sub_.append(mul(qm, Pq, p))
        fams['C[Q,K,H,diag]'] = sub_
        # (Q)^k * (forms of degree d-2k, even)
        for k in (d//2 - 1, d//2 - 2):
            if k <= 0: continue
            Qm = poly_pow_monos(Qp, k, p)
            lst = []
            for qm in Qm.values():
                for m in monomials(d - 2*k, 0):
                    lst.append(mul(qm, {m: 1}, p))
            fams[f'(Q)^{k}*S_{d-2*k}'] = lst

    if par == 1 and d % 2 == 0 and d >= 8:
        oddq = [{m: 1} for m in monomials(2, 1)]  # X_i Y_j
        Qm = poly_pow_monos(Qp, (d-2)//2, p)
        fams['oddquad*G_k(Q)'] = [mul(o, qm, p) for o in oddq for qm in Qm.values()]
        quartics = Kp + Hp + diag
        A0 = reduce_mod(P('X0*Y0'), p)
        lst = []
        for nq in range(0, (d-2)//4 + 1):
            dq = d - 2 - 4*nq
            if dq < 0 or dq % 2: continue
            Qm2 = poly_pow_monos(Qp, dq//2, p)
            for combo in itertools.combinations_with_replacement(range(9), nq):
                Pq = A0
                for c in combo: Pq = mul(Pq, quartics[c], p)
                for qm in Qm2.values(): lst.append(mul(qm, Pq, p))
        fams['X0Y0*C[Q,K,H,diag]'] = lst
    # plane pair vanishing order k = d//2 (products of cones over Pi and Pi-bar): symmetric (even) / antisym (odd)
    if d % 2 == 0:
        k = d // 2
        u, ub = plane_forms((1,1,1), p)
        Gu = poly_pow_monos(u, k, p); Gub = poly_pow_monos(ub, k, p)
        keys = sorted(Gu)
        lst = []
        for a in range(len(keys)):
            for b in range(a, len(keys)):
                t1 = mul(Gu[keys[a]], Gub[keys[b]], p); t2 = mul(Gu[keys[b]], Gub[keys[a]], p)
                s = add(t1, t2, p) if par == 0 else sub(t1, t2, p)
                if s: lst.append(s)
        fams[f'I_Pi^{k} cap I_Pibar^{k}'] = lst
    # keep only polys of the right parity / degree
    out = {}
    for name, lst in fams.items():
        good = [P_ for P_ in lst if P_ and core.degree(P_) == d and parity_of(P_) == par]
        out[name] = good
    return out

def J_basis(d, par, p):
    Fp = reduce_mod(F, p); Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    lst = []
    for i in range(6):
        parM = (par - (0 if i < 3 else 1)) % 2
        for m in monomials(d - 5, parM):
            lst.append(core.shift(Fi[i], m, 1, p))
    return lst

def rank_of(polys, p, monos=None):
    if not polys: return 0, None
    if monos is None:
        monos = {}
        for P_ in polys:
            for e in P_: monos.setdefault(e, len(monos))
    R, C = len(monos), len(polys)
    flat = [0]*(R*C)
    for j, P_ in enumerate(polys):
        for e, v in P_.items():
            flat[monos[e]*C + j] = v % p
    return nmod_mat(R, C, flat, p).rank(), monos

def intersect_dims(kernel, fam, J, p):
    """dims: fam, fam∩ker, fam∩J, (fam∩ker) mod J  (all via ranks)"""
    monos = {}
    for lst in (kernel, fam, J):
        for P_ in lst:
            for e in P_: monos.setdefault(e, len(monos))
    rk = lambda L: rank_of(L, p, monos)[0] if L else 0
    rK, rF, rJ = rk(kernel), rk(fam), rk(J)
    rKF = rk(kernel + fam); rJF = rk(J + fam); rKJ = rk(kernel + J); rKJF = rk(kernel + J + fam)
    dimFK = rF + rK - rKF
    dimFJ = rF + rJ - rJF
    # (F ∩ K) mod J : dim = dim(F∩K) - dim(F∩K∩J); F∩K∩J = F∩J (since J ⊂ K)
    return {'dim_fam': rF, 'dim_ker': rK, 'dim_J': rJ, 'dim_ker_plus_J': rKJ, 'fam_cap_ker': dimFK, 'fam_cap_J': dimFJ, 'fam_cap_ker_mod_J': dimFK - dimFJ}

if __name__ == '__main__':
    d = int(sys.argv[1]); par = int(sys.argv[2]); p = core.PRIMES[1] if len(sys.argv) < 4 else int(sys.argv[3])
    import pickle, os
    pk = f'{WD}/e1kernel_basis_d{d}_par{par}.pkl'
    if os.path.exists(pk):
        kernel = pickle.load(open(pk, 'rb'))
    else:
        # text basis (coeffs centred mod PRIMES[0]) -- only valid if p == PRIMES[0]
        assert p == core.PRIMES[0]
        kernel = [core.parse(l.strip(), p) for l in open(f'{WD}/e1kernel_basis_d{d}_par{par}.txt') if l.strip()]
        kernel = [reduce_mod(k_, p) for k_ in kernel]
    J = J_basis(d, par, p)
    fams = families(d, par, p)
    for name, fam in fams.items():
        info = intersect_dims(kernel, fam, J, p)
        print(d, par, name, json.dumps(info), flush=True)
