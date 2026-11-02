"""Exact rational basis of the primitive E1 kernel  Ann([A])_d  in R_d = S_d/J(F)_d  (one parity), by multi-prime
CRT + rational reconstruction of the canonical (rref) basis, followed by EXACT verification over Q that
A * D_k = sum_i U_{k,i} F_i  (cofactors lifted the same way).  Receipt: kernel_d<d>_par<par>_exact.json.
usage: python3 kernel_lift.py d parity [maxprimes]"""
import sys, time, json
from fractions import Fraction as Fr
from flint import nmod_mat
from common import *
from hodge import jacobian
from hodge.core import monomials, shift, pderiv, mul, reduce_mod
F, A, q, rec = core.load_instance(IID)
d = int(sys.argv[1]); par = int(sys.argv[2]); maxp = int(sys.argv[3]) if len(sys.argv) > 3 else 10
Fi = [pderiv(F, i) for i in range(6)]
t0 = time.time()
Dm = monomials(d, par); K = len(Dm)
tpar = (1 + par) % 2
gens = [(i, m) for i in range(6) for m in monomials(d + 1, (tpar - (0 if i < 3 else 1)) % 2)]
rows = monomials(d + 6, tpar); ridx = {e: k for k, e in enumerate(rows)}
R, C = len(rows), K + len(gens)
print(f'd={d} par={par}: matrix {R} x {C} (D cols {K}, J gens {len(gens)})', flush=True)

def kernel_basis_modp(p):
    """Canonical rref basis (mod p) of {D in S_d^par : A D in J}, projected to D-coordinates, then reduced to a
    canonical basis of the quotient by J_d (standard monomials), returned with the pivot structure."""
    Ap = reduce_mod(A, p); Fip = [reduce_mod(f, p) for f in Fi]
    flat = [0] * (R * C)
    for j, m in enumerate(Dm):
        for e, v in shift(Ap, m, 1, p).items(): flat[ridx[e] * C + j] = v
    for j, (i, m) in enumerate(gens):
        for e, v in shift(Fip[i], m, 1, p).items(): flat[ridx[e] * C + K + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], p)
    Pr, rk = proj.rref()
    vecs = [{Dm[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])} for i in range(rk)]
    red = jacobian._Reducer(reduce_mod(F, p), d, p, par)
    sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
    nf = [red.normal_form(v) for v in vecs]
    flat2 = [0] * (len(nf) * len(sm))
    for r_, v in enumerate(nf):
        for e, c in v.items(): flat2[r_ * len(sm) + sidx[e]] = c
    Q = nmod_mat(len(nf), len(sm), flat2, p)
    Qr, rk2 = Q.rref()
    basis = [{sm[j]: int(Qr[i, j]) for j in range(len(sm)) if int(Qr[i, j])} for i in range(rk2)]
    piv = []
    for i in range(rk2):
        j = 0
        while int(Qr[i, j]) == 0: j += 1
        piv.append(sm[j])
    return basis, piv, rk, red.jacobian_dim

cur = None; Mprod = 1; used = []
for p in core.PRIMES[:maxp]:
    basis, piv, kdim, jdim = kernel_basis_modp(p)
    print(f'p={p}: kernel {kdim}, J_d {jdim}, primitive {len(basis)}, pivots {piv[:3]}...', round(time.time() - t0, 1), flush=True)
    if cur is None:
        cur = basis; pivref = piv; Mprod = p; used = [p]; nb = len(basis); keys = [set(b) for b in basis]
    else:
        if len(basis) != nb or piv != pivref:
            print('pivot structure differs mod', p, '- skipping prime', flush=True); continue
        new = []
        for b_old, b_new in zip(cur, basis):
            allk = set(b_old) | set(b_new)
            new.append({e: core.crt_pair(b_old.get(e, 0), Mprod, b_new.get(e, 0), p) for e in allk})
        cur = new; Mprod *= p; used.append(p)
    # try reconstruction
    try:
        rat = [{e: core.ratrecon(v, Mprod) for e, v in b.items()} for b in cur]
    except Exception as ex:
        print('ratrecon exception', ex, flush=True); continue
    if any(v is None for b in rat for v in b.values()):
        print('ratrecon incomplete with', len(used), 'primes', flush=True); continue
    rat = [{e: Fr(v) for e, v in b.items() if v} for b in rat]
    # exact verification: A*D_k in J over Q via rational cofactors (solve mod primes with fixed pivots, lift, verify)
    ok_all = True; cof_all = []
    for k_, Dk in enumerate(rat):
        AD = mul(A, Dk)
        # cofactor system: rows = monomials deg d+6 parity tpar; cols = gens
        Cg = len(gens)
        sols = []; mods = []; piv_cols = None; rpiv = None
        ok = False
        for p2 in core.PRIMES[:maxp]:
            Fip = [reduce_mod(f, p2) for f in Fi]
            flat = [0] * (R * Cg)
            for j, (i, m) in enumerate(gens):
                for e, v in shift(Fip[i], m, 1, p2).items(): flat[ridx[e] * Cg + j] = v
            G = nmod_mat(R, Cg, flat, p2)
            b = [0] * R
            for e, v in reduce_mod(AD, p2).items(): b[ridx[e]] = v
            if piv_cols is None:
                Gr, rkG = G.rref()
                piv_cols = []; jj = 0
                for i_ in range(rkG):
                    while jj < Cg and int(Gr[i_, jj]) == 0: jj += 1
                    piv_cols.append(jj); jj += 1
                GT = nmod_mat(Cg, R, [int(G[r_, c_]) for c_ in range(Cg) for r_ in range(R)], p2)  # transpose
                sub = nmod_mat(R, rkG, [int(G[r_, c_]) for r_ in range(R) for c_ in piv_cols], p2)
                ST, rkT = sub.transpose().rref()
                rpiv = []; jj = 0
                for i_ in range(rkT):
                    while jj < R and int(ST[i_, jj]) == 0: jj += 1
                    rpiv.append(jj); jj += 1
            Sq = nmod_mat(len(rpiv), len(piv_cols), [int(G[r_, c_]) for r_ in rpiv for c_ in piv_cols], p2)
            try:
                x = Sq.solve(nmod_mat(len(rpiv), 1, [b[r_] for r_ in rpiv], p2))
            except Exception:
                continue
            xs = [int(x[i_, 0]) for i_ in range(len(piv_cols))]
            sols.append(xs); mods.append(p2)
            Mp = 1; cc = None
            for xs_, pp in zip(sols, mods):
                if cc is None: cc = list(xs_); Mp = pp
                else: cc = [core.crt_pair(a_, Mp, c_, pp) for a_, c_ in zip(cc, xs_)]; Mp *= pp
            rr = [core.ratrecon(a_, Mp) for a_ in cc]
            if any(v is None for v in rr): continue
            U = [dict() for _ in range(6)]
            for t_, c_ in enumerate(piv_cols):
                if rr[t_]:
                    i, m = gens[c_]; U[i][m] = Fr(rr[t_])
            tot = {}
            for i in range(6):
                for e, v in mul(U[i], Fi[i]).items(): tot[e] = tot.get(e, 0) + v
            if all(Fr(tot.get(e, 0)) == Fr(AD.get(e, 0)) for e in set(tot) | set(AD)):
                ok = True; cof_all.append([core.to_str(u) for u in U]); break
        print(f'  basis element {k_}: exact A*D in J over Q: {ok} ({len(mods)} primes)', round(time.time() - t0, 1), flush=True)
        if not ok: ok_all = False; break
    if ok_all:
        receipt(f'kernel_d{d}_par{par}_exact.json', {'d': d, 'parity': par, 'dim': len(rat), 'primes_used': used,
                 'statement': 'canonical rref basis of Ann([A])_d (standard monomials of R_d), each verified A*D = sum U_i F_i over Q',
                 'basis': [core.to_str(b) for b in rat], 'cofactors': cof_all})
        log({'exp': 'e1-kernel-exact-lift', 'd': d, 'parity': par, 'dim': len(rat), 'primes_used': used,
             'all_verified_over_Q': True, 'seconds': round(time.time() - t0, 1)})
        break
else:
    log({'exp': 'e1-kernel-exact-lift', 'd': d, 'parity': par, 'all_verified_over_Q': False, 'note': 'lift failed', 'seconds': round(time.time() - t0, 1)})
