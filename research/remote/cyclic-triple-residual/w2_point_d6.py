"""Identify the (at most one) nonspurious d=6 Y-even relaxation point: numerical relaxation search, rational
reconstruction, exact relaxation membership (Bil(D, N0) in J_12 over Q), exact full-E2 test. usage: w2_point_d6.py [restarts]"""
import os, sys, time, json
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from w2core import *
from flint import fmpz_mat
import numpy as np
n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
F, A, q = load(); G = full_group(F, A); S = Sect(F, A, G)
Sp = Sect(F, A, parity_group())
r = e1_space(Sp, F, A, q, 6, (0,), verbose=False)
fam = Family(Sp, F, A, q, 6, (0,), r['basis'], verbose=False)
t0 = time.time()
cands = fam.search(n_restarts=n, seed=23, relax=True, maxiter=500)
rels = sorted(c['rel'] for c in cands)
print('relax rels (sorted, first 10):', [f'{x:.2e}' for x in rels[:10]], f'({time.time()-t0:.1f}s)', flush=True)
best = cands[0]
Df = fam.D_float(best['a'])
print('best D (float, |coef| sorted):', sorted(((round(v, 8), core.to_str({m: 1})) for m, v in Df.items()), key=lambda t: -abs(t[0]))[:20], flush=True)
Fz = sect.integral(F); dFz = [sect.ideriv(Fz, i) for i in range(6)]
Fi = [pderiv(F, i) for i in range(6)]
def exact_relax(D):
    """Is Bil(D, N0) in J_12 over Q, for D = sum w_i F_i (w recovered exactly by solving the E1 system)? Use any exact lift."""
    # exact lift via the E1 nullspace route: solve sum N_i F_i = A D over Q in the G-sector (small)
    Dz = sect.integral(D); c = sect.chi_of(Dz, S.chi)
    # columns m' F_i (deg 7) in sector chiA + c; rhs A D
    cols = []; meta = []
    for i in range(6):
        for mp in S.mons(7, S.chiN(c, i)):
            cols.append(sect.imul({mp: 1}, dFz[i])); meta.append((i, mp))
    rhs = sect.imul(sect.integral(A), Dz)
    rows = {}
    for P in cols + [rhs]:
        for e in P: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, P in enumerate(cols):
        for e, v in P.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    from flint import fmpq_mat
    M = fmpq_mat(fmpz_mat(R, C + 1, flat))
    Mr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = 0
        while j <= C and Mr[i, j] == 0: j += 1
        piv.append(j)
    if piv and piv[-1] == C:
        return None, 'E1 fails exactly'
    N = [{} for _ in range(6)]
    LF = Fraction(Fz[next(iter(Fz))]) / Fraction(F[next(iter(Fz))])
    for i_, j in enumerate(piv):
        v = Mr[i_, C]; v = Fraction(int(v.p), int(v.q))
        if v:
            i, mp = meta[j]
            N[i][mp] = v * LF / q    # columns used Fz_i = LF * F_i ; A D (integral A, Dz) ... rescale below
    # rescale: we solved sum N_i (LF F_i) = A Dz  with Dz = LD * D  ->  sum (N_i LF / LD) F_i = A D
    LD = Fraction(Dz[next(iter(Dz))]) / Fraction(D[next(iter(Dz))])
    N = [{e: v / LD for e, v in Ni.items()} for Ni in N]
    assert not e1_check(F, A, q, D, N), 'lift check failed'
    B = bil(D, N)
    if not B:
        return True, 'Bil = 0'
    cB = sect.chi_of(sect.integral(B), S.chi)
    gens = [sect.imul({m: 1}, dFz[i]) for i in range(6) for m in S.mons(7, S.chiJgen(cB, i))]
    Bz = sect.integral(B)
    rows2 = {}
    for g in gens + [Bz]:
        for e in g: rows2.setdefault(e, len(rows2))
    def mat(polys):
        R2, C2 = len(rows2), len(polys); fl = [0] * (R2 * C2)
        for j, P in enumerate(polys):
            for e, v in P.items(): fl[rows2[e] * C2 + j] = v
        return fmpz_mat(R2, C2, fl)
    r0 = mat(gens).rank(); r1 = mat(gens + [Bz]).rank()
    return r0 == r1, f'rank J {r0}, with B {r1}'
out = []
for D, descr, err in fam.reconstruct(best['a'], dens=(1, 2, 3, 4, 6, 7, 12, 14, 21, 28, 42, 49, 84, 98, 147, 196, 294, 588, 1176, 2401, 4802, 10**4, 10**5, 10**6)):
    if core.degree(D) != 6: continue
    rec = {'descr': descr, 'recon_err': err, 'terms': len(D), 'D': core.to_str(D)[:400]}
    try:
        ok, info = exact_relax(D)
        rec['relax_exact'] = ok; rec['relax_info'] = info
        if ok:
            conf = confirm_exact(F, A, q, D, exact=True)
            rec['fullE2'] = {k: conf.get(k) for k in ('e1', 'consistent', 'rank', 'rank_aug', 'cols', 'rows', 'verified')}
            if conf.get('witness'):
                rec['WITNESS'] = True
                json.dump({'instanceId': IID, **conf['witness']}, open(f'{WD}/witness-candidate.json', 'w'), indent=1)
    except Exception as ex:
        rec['error'] = repr(ex)
    print(rec, flush=True)
    out.append(rec)
    if rec.get('relax_exact'):
        break
w2log({'batch': 'relaxation-point-d6-identify', 'd': 6, 'parity': 0, 'restarts': n, 'relax_rels_sorted_head': rels[:10],
       'best_D_float': {core.to_str({m: 1}): v for m, v in Df.items()}, 'reconstructions': out, 'seconds': round(time.time() - t0, 1)})
print('done', round(time.time() - t0, 1), 's', flush=True)
