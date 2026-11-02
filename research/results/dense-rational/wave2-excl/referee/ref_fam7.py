#!/usr/bin/env python3
"""Referee's own generator: degree-7 mixed-parity Hessian relaxation restricted to a structured linear family of polar
fields w (D = w(F), w quadratic, any Y-parity), condition w^T Hess(F) w in J_8 (both parities), gauge w -> w + l*Euler
(l linear) removed inside the family. Families: weuler (w_i = z_i g_i, g_i linear), blockdiag (w_X in C[X]_2, w_Y in C[Y]_2),
cross (w_X in C[Y]_2, w_Y in C[X]_2), gradient (w = grad g, g cubic). The family and gauge vectors are integer vectors; the
gauge-in-family intersection is computed both over Q (exact rank formula) and mod p and must agree.
Usage: ref_fam7.py prime fam1,fam2,..."""
import subprocess
from ref_common import *
from ref_ltparse import check
p = int(sys.argv[1]); fams = sys.argv[2].split(',')
F, A = load(); Fp = modp(F, p)
Fi = [pd(Fp, i, p) for i in range(NV)]
H = [[pd(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
t0 = time.time()
mons8 = monos(8); i8 = {m: j for j, m in enumerate(mons8)}
gens = [shift(Fi[i], m, p) for i in range(NV) for m in monos(3)]
G = nmod_mat(len(gens), len(mons8), p); fill_rows(G, gens, i8)
Gr, rk = G.rref(); assert rk == 336 == len(gens), rk
piv = pivots_of(Gr, rk); pivset = set(piv); nonpiv = [j for j in range(len(mons8)) if j not in pivset]
assert len(nonpiv) == 951
basis2 = monos(2); b2 = {m: k for k, m in enumerate(basis2)}            # coordinates: (i, m) -> 21*i + b2[m]
def coords(field):
    v = [0] * (NV * len(basis2))
    for i in range(NV):
        for m, c in field[i].items(): v[NV * 0 + len(basis2) * i + b2[m]] = c % p
    return v
def unit(i, m): f = [{} for _ in range(NV)]; f[i] = {m: 1}; return f
ei = lambda i: tuple(1 if j == i else 0 for j in range(NV))
gauge = []
for l in monos(1):
    f = [{tuple(a + b for a, b in zip(l, ei(i))): 1} for i in range(NV)]; gauge.append(f)
def family(fam):
    out = []
    if fam == 'weuler':
        for i in range(NV):
            for g in monos(1): out.append(unit(i, tuple(a + b for a, b in zip(g, ei(i)))))
    elif fam == 'blockdiag':
        for i in range(NV):
            for m in basis2:
                if (m[3] + m[4] + m[5] == 0) if i < 3 else (m[0] + m[1] + m[2] == 0): out.append(unit(i, m))
    elif fam == 'cross':
        for i in range(NV):
            for m in basis2:
                if (m[0] + m[1] + m[2] == 0) if i < 3 else (m[3] + m[4] + m[5] == 0): out.append(unit(i, m))
    elif fam == 'gradient':
        for g in monos(3): out.append([pd({g: 1}, i, p) for i in range(NV)])
    else: raise SystemExit('unknown family ' + fam)
    return out
print(f'J_8 quotient 951 ready ({time.time()-t0:.1f}s)', flush=True)
for fam in fams:
    t0 = time.time(); fields = family(fam); nf = len(fields)
    Vf = [coords(f) for f in fields]; Gv = [coords(f) for f in gauge]
    n = NV * len(basis2)
    Mf = nmod_mat(n, nf, [Vf[c][r] for r in range(n) for c in range(nf)], p); assert Mf.rank() == nf
    # gauge-in-family: exact over Q via integer ranks (entries are small integers, fmpq exact), and mod p via nullspace
    def qmat(cols):
        M = fmpq_mat(n, len(cols))
        for c, v in enumerate(cols):
            for r in range(n):
                if v[r]: M[r, c] = fmpq(v[r] if v[r] <= p // 2 else v[r] - p)
        return M
    rF, rG, rFG = qmat(Vf).rank(), qmat(Gv).rank(), qmat(Vf + Gv).rank()
    dim_int_Q = rF + rG - rFG
    MK = nmod_mat(n, nf + len(Gv), [Vf[c][r] if c < nf else (p - Gv[c - nf][r]) % p for r in range(n) for c in range(nf + len(Gv))], p)
    K, nk = MK.nullspace(); assert nk == dim_int_Q, (nk, dim_int_Q)
    gin = nmod_mat(max(nk, 1), nf, [int(K[j, k]) for k in range(nk) for j in range(nf)] if nk else [0] * nf, p)
    ginr, grk = gin.rref(); drop = set(pivots_of(ginr, grk)) if nk else set()
    keep = [a for a in range(nf) if a not in drop]
    # pair products P_ab = w_a^T Hess w_b, normal form mod J_8
    pairs = [(a, b) for a in range(nf) for b in range(a, nf)]
    V = nmod_mat(len(pairs), len(mons8), p)
    for r, (a, b) in enumerate(pairs):
        P = {}
        for i in range(NV):
            if not fields[a][i]: continue
            for j in range(NV):
                if not fields[b][j]: continue
                P = add(P, mul(mul(fields[a][i], fields[b][j], p), H[i][j], p), p)
        for mm, c in P.items(): V[r, i8[mm]] = c
    VP = nmod_mat(len(pairs), rk, [int(V[r, j]) for r in range(len(pairs)) for j in piv], p)
    NF = V - VP * Gr
    eqs = []
    for j in nonpiv:
        terms = []
        for r, (a, b) in enumerate(pairs):
            if a in drop or b in drop: continue
            c = int(NF[r, j])
            if a != b: c = 2 * c % p
            if c: terms.append(f'{c}*c{a}^2' if a == b else f'{c}*c{a}*c{b}')
        if terms: eqs.append('+'.join(terms))
    fn = os.path.join(HERE, f'ref_fam7_mixed_{fam}_p{p}.ms'); out = fn.replace('.ms', '.g1.out')
    with open(fn, 'w') as f: f.write(','.join(f'c{a}' for a in keep) + '\n' + f'{p}\n' + ',\n'.join(eqs) + '\n')
    tb = time.time() - t0
    r = subprocess.run(['msolve', '-v', '0', '-t', '4', '-g', '1', '-f', fn, '-o', out], capture_output=True, text=True, timeout=3600)
    rec = check(out)
    log(check='referee-fam7-own-generator', family=fam, sector='mixed', d=7, prime=p, family_dim=nf, gauge_in_family_Q=dim_int_Q,
        gauge_in_family_modp=nk, unknowns=len(keep), equations=len(eqs), build_seconds=round(tb, 1), msolve_exit=r.returncode,
        **{k: rec[k] for k in rec if k not in ('file', 'prime', 'nvars')},
        verdict=('EMPTY mod p (m-primary) -> no legal degree-7 D = w(F) with w in this family, over C' if rec.get('m_primary_mod_p') else 'NOT CERTIFIED'),
        seconds=round(time.time() - t0, 1))
