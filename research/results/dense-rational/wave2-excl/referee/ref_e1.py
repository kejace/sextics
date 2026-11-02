#!/usr/bin/env python3
"""Referee's own E1 bookkeeping at a fresh prime, with BOTH rank caps made explicit:
  K_full(d,parD) = {D in S_d^{parD} : A D in J_{d+6}}  via nullspace of [A*S_d | S_{d+1}*grad F] (rows S_{d+6}^{1-parD});
  rigorous over Q iff rank_p(N-block) = dim J_{d+6} (exact Hilbert count) AND rank_p(J_d-block) = dim J_d AND dim K_full = dim J_d.
Also: exact fmpq ranks (A not in J_6^-, Sym^2 -> S_4 injective) and a mod-p check of the relaxation identity.
Usage: ref_e1.py prime d1,d2,...  [--skip-exact]"""
import random
from ref_common import *
p = int(sys.argv[1]); ds = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 and sys.argv[2] else []
F, A = load(); S, R = hilbert()
FiQ = [pd(F, i) for i in range(NV)]
if '--skip-exact' not in sys.argv:
    t0 = time.time()
    m6 = monos(6, 1); i6 = {m: j for j, m in enumerate(m6)}
    gens = [shift(FiQ[i], m) for i in range(NV) for m in monos(1, (1 - pi(i)) % 2)]
    def fq(P, M, j):
        for e, c in P.items():
            fr = Fraction(c); M[i6[e], j] = fmpq(fr.numerator, fr.denominator)
    M0 = fmpq_mat(len(m6), len(gens)); M1 = fmpq_mat(len(m6), len(gens) + 1)
    for j, g in enumerate(gens): fq(g, M0, j); fq(g, M1, j)
    fq(A, M1, len(gens))
    log(check='referee-exact-A-not-in-J6', field='Q exact (fmpq_mat)', rank_J6_minus=M0.rank(), rank_aug=M1.rank(), n_gens=len(gens),
        verdict=('A NOT in J_6 over Q' if M1.rank() > M0.rank() else 'A in J_6'), seconds=round(time.time() - t0, 2))
    t0 = time.time()
    Fij = [[pd(FiQ[i], j) for j in range(NV)] for i in range(NV)]
    m4 = monos(4); i4 = {m: j for j, m in enumerate(m4)}
    prs = [(i, j) for i in range(NV) for j in range(i, NV)]
    M = fmpq_mat(len(m4), len(prs))
    for c, (i, j) in enumerate(prs):
        for e, v in Fij[i][j].items():
            fr = Fraction(v) * (1 if i == j else 2); M[i4[e], c] = fmpq(fr.numerator, fr.denominator)
    log(check='referee-exact-sym2-to-S4-rank', field='Q exact (fmpq_mat)', rows=len(m4), cols=len(prs), rank=M.rank(),
        verdict=('injective: v^T Hess(F) v = 0 forces v = 0' if M.rank() == len(prs) else 'KERNEL'), seconds=round(time.time() - t0, 2))
    # relaxation identity mod p for a random mixed-parity LINEAR w (D = w(F) of degree 6, N = A w)
    t0 = time.time(); random.seed(11)
    Fp = modp(F, p); Ap = modp(A, p); Fi = [pd(Fp, i, p) for i in range(NV)]
    H = [[pd(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    w = [{m: random.randrange(1, p) for m in monos(1)} for i in range(NV)]
    D = {}
    for i in range(NV): D = add(D, mul(w[i], Fi[i], p), p)
    N = [mul(Ap, w[i], p) for i in range(NV)]
    e1 = {}
    for i in range(NV): e1 = add(e1, mul(N[i], Fi[i], p), p)
    assert not sub(e1, mul(Ap, D, p), p)
    divN = {}
    for i in range(NV): divN = add(divN, pd(N[i], i, p), p)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N[i], pd(D, i, p), p), p)
    B = sub(mul(D, divN, p), ND, p)
    wHw = {}
    for i in range(NV):
        for j in range(NV): wHw = add(wHw, mul(mul(w[i], w[j], p), H[i][j], p), p)
    AwHw = mul(Ap, wHw, p)
    m12 = monos(12); i12 = {m: j for j, m in enumerate(m12)}
    gens12 = [shift(Fi[i], m, p) for i in range(NV) for m in monos(7)]
    def in_J12(v):
        G0 = nmod_mat(len(m12), len(gens12), p); fill(G0, gens12, i12)
        G1 = nmod_mat(len(m12), len(gens12) + 1, p); fill(G1, gens12 + [v], i12)
        r0 = G0.rank(); r1 = G1.rank(); return r0, r1
    r0, rp = in_J12(add(B, AwHw, p)); _, rm = in_J12(sub(B, AwHw, p))
    log(check='referee-relaxation-identity-modp', prime=p, w='random mixed-parity linear field (seed 11)', rank_J12=r0,
        rank_with_B_plus_AwHw=rp, rank_with_B_minus_AwHw=rm,
        verdict=('B_D(Aw) == -A w^T Hess(F) w  (mod J)' if (rp == r0 and rm > r0) else 'IDENTITY NOT CONFIRMED'), seconds=round(time.time() - t0, 1))
Fp = modp(F, p); Ap = modp(A, p); Fi = [pd(Fp, i, p) for i in range(NV)]
for d in ds:
    for parD in (0, 1):
        t0 = time.time(); parT = 1 - parD
        Dm = monos(d, parD); nD = len(Dm)
        rows = monos(d + 6, parT); ridx = {e: j for j, e in enumerate(rows)}
        colsA = [shift(Ap, m, p) for m in Dm]
        colsN = [shift(Fi[i], m, p) for i in range(NV) for m in monos(d + 1, (parT - pi(i)) % 2)]
        Rn, Cn = len(rows), nD + len(colsN)
        M = nmod_mat(Rn, Cn, p); fill(M, colsA + colsN, ridx)
        MN = nmod_mat(Rn, len(colsN), p); fill(MN, colsN, ridx)
        rankN = MN.rank(); del MN
        X, nullity = M.nullspace(); rank = Cn - nullity
        proj = nmod_mat(nullity, nD, [int(X[j, k]) for k in range(nullity) for j in range(nD)], p) if nullity else None
        dimK = proj.rank() if nullity else 0
        Jg = [shift(Fi[i], m, p) for i in range(NV) for m in monos(d - 5, (parD - pi(i)) % 2)]
        Di = {m: j for j, m in enumerate(Dm)}
        rJ = 0
        if Jg:
            JM = nmod_mat(len(Jg), nD, p); fill_rows(JM, Jg, Di); rJ = JM.rank()
        exJd = S[d][parD] - R[d][parD]; exJd6 = S[d + 6][parT] - R[d + 6][parT]
        forced = max(0, R[d][parD] - R[d + 6][parT])
        rig = (rankN == exJd6) and (rJ == exJd) and (dimK - exJd == forced)
        log(check='referee-e1-kernel', d=d, parD=parD, prime=p, rows=Rn, cols=Cn, nD=nD, rank_full=rank, rank_Nblock=rankN,
            exact_dim_J_d6=exJd6, Ncap_attained=(rankN == exJd6), dim_Kfull_modp=dimK, rank_Jd_block=rJ, exact_dim_Jd=exJd,
            Jd_cap_attained=(rJ == exJd), kernel_mod_J=dimK - exJd, forced_lower_bound=forced, rigorous_over_Q=rig,
            verdict=(f'(J:A)_{d} = J_{d} over Q (both caps attained)' if rig and dimK == exJd else
                     (f'kernel mod J = {dimK - exJd}, exact over Q' if rig else 'MODULAR ONLY (a cap not attained)')),
            seconds=round(time.time() - t0, 1))
