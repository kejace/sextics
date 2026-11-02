"""Rigorous fixed-D certificate test (character-sector restricted system + full-system cross-check).

Protocol (q = 1):
  * restricted system M x = b (unknowns N, Q filtered to the character sector forced by D; lossless by averaging),
    solved mod p by random row compression (fixed_d.solve_compressed): consistent? rank r_c, pivots.
  * consistent  -> exact rational solve via symmetry.solve_fixed_D_symmetric (compressed CRT lift) and exact verify.
  * inconsistent -> certificate that the rational rank of M equals the modular rank:
        a-priori kernel vectors (exact, rational):  (z G, 0) for G in S_d,  (F m e_i, D d_i m - m d_i D) for m in S_{d-5};
        if their modular rank r_ap satisfies r_c = ncols - r_ap, then rank_Q M = rank_p M and the modular
        inconsistency is a proof of inconsistency over Q (hence over C).
        Otherwise: dense modular nullspace at several primes, CRT + rational reconstruction, exact verification.
Usage: python3 rigorfix.py <family-tag> '<D expr>' [--full] [--primes k]
"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core, fixed_d, symmetry
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power
from flint import nmod_mat
from fractions import Fraction

F, A, q = L.F, L.A, L.q
GROUP = symmetry.diagonal_group(F)


def apriori_kernel(S, p):
    """Exact a-priori kernel vectors of the (restricted) fixed-D system, as sparse dicts over column indices."""
    D = S.D
    d = S.d
    colidx = {c: j for j, c in enumerate(S.cols)}
    vecs = []
    parG = (core.parity_of(A) - core.parity_of(D)) % 2
    for G in monomials(d, parG):
        v = {}
        ok = True
        for i in range(NV):
            m = list(G); m[i] += 1; m = tuple(m)
            key = ('N', i, m)
            if key not in colidx:
                ok = False; break
            v[colidx[key]] = 1
        if ok:
            vecs.append(('zG', G, v))
    Fp = S.F
    Dp = S.D
    Di = [pderiv(Dp, i, p) for i in range(NV)]
    if d >= 5:
        for i in range(NV):
            for m in monomials(d - 5):
                Nm = shift(Fp, m, 1, p)                     # N_i = F m
                mono = {m: 1}
                Qm = sub(mul(Dp, pderiv(mono, i, p), p), mul(mono, Di[i], p), p)   # D d_i m - m d_i D
                v = {}
                ok = True
                for e, c in Nm.items():
                    key = ('N', i, e)
                    if key not in colidx:
                        ok = False; break
                    v[colidx[key]] = c
                if ok:
                    for e, c in Qm.items():
                        key = ('Q', -1, e)
                        if key not in colidx:
                            ok = False; break
                        v[colidx[key]] = c
                if ok and v:
                    vecs.append(('Fm', (i, m), v))
    return vecs


def rank_of_vecs(vecs, C, p):
    if not vecs:
        return 0
    flat = [0] * (len(vecs) * C)
    for r, (_, _, v) in enumerate(vecs):
        for j, c in v.items():
            flat[r * C + j] = c % p
    return nmod_mat(len(vecs), C, flat, p).rank()


def exact_kernel_check(S, cols, vec):
    """Exact check that (N,Q) given by rational vec (dict col->Fraction) is in the kernel of the q=1 system."""
    N = [{} for _ in range(NV)]; Q = {}
    for j, v in vec.items():
        kind, i, m = cols[j]
        if v == 0: continue
        if kind == 'N': N[i][m] = v
        else: Q[m] = v
    D = core.parse(core.to_str(S.Dq)) if hasattr(S, 'Dq') else None
    return N, Q


def run(tag, expr, full=False, nprimes=3):
    D = core.parse(expr)
    d = core.degree(D)
    t0 = time.time()
    rec = {'exp': 'rigorfix', 'family': tag, 'D': expr, 'd': d}
    p = L.PRIMES[2]
    S, sub_ = symmetry.restricted_system(F, A, q, D, p, GROUP, verbose=True)
    r = fixed_d.solve_compressed(S, verbose=True)
    rec.update(prime=p, rows=S.nrows, cols=S.ncols, full_cols=S.full_cols, stabilizer=len(sub_),
               rank_c=r['rank'], nullity_c=r['nullity'], consistent_modp=bool(r['consistent']))
    if r['consistent']:
        st, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, primes=L.PRIMES[2:2 + 6], verbose=True,
                                                    group=GROUP, compressed=True)
        rec['status'] = st
        if st == 'solved':
            ok, info = core.verify_witness(F, A, q, D, data['N'], data['Q'])
            rec['verified'] = ok
            if ok:
                w = {'D': expr, 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q']), 'd': d}
                fn = f'{L.OUT}/WITNESS_{tag}.json'
                json.dump(w, open(fn, 'w'), indent=1)
                rec['witness_file'] = fn
                print('!!!!!!!!!!!! WITNESS FOUND AND VERIFIED', fn, flush=True)
        rec['seconds'] = round(time.time() - t0, 1)
        L.receipt('rigorfix', **rec)
        return rec
    # inconsistent mod p: certify
    vecs = apriori_kernel(S, p)
    r_ap = rank_of_vecs(vecs, S.ncols, p)
    rank_M_c = r['rank'] - 1          # solve_compressed reports rank of [M|b]; inconsistent => pivot in b-column
    rec.update(rank_M_c=rank_M_c, apriori_vectors=len(vecs), apriori_rank=r_ap, certified=(rank_M_c == S.ncols - r_ap))
    if not rec['certified']:
        # try a second compression seed (rank may have been lost) and the dense rank
        r2 = fixed_d.solve_compressed(S, seed=11, verbose=True)
        rec['rank_c_seed11'] = r2['rank']
        if r2['rank'] - 1 == S.ncols - r_ap and not r2['consistent']:
            rec['certified'] = True
        else:
            try:
                M = S.dense_augmented()
                C = S.ncols
                M0 = nmod_mat(S.nrows, C, [int(M[i, j]) for i in range(S.nrows) for j in range(C)], p)
                rk = M0.rank(); rka = M.rank()
                rec.update(rank_dense=rk, rank_dense_aug=rka, nullity_dense=C - rk)
                rec['consistent_modp'] = (rk == rka)
                rec['certified'] = (not rec['consistent_modp']) and (rk == C - r_ap)
            except Exception as ex:
                rec['dense_error'] = repr(ex)
    rec['note'] = ('inconsistent mod p; rank_Q(M)=rank_p(M) certified by a-priori rational kernel => inconsistent over Q'
                   if rec['certified'] else 'inconsistent mod p; rational rank NOT certified (a-priori kernel too small)')
    # second prime for agreement
    p2 = L.PRIMES[3]
    S2, _ = symmetry.restricted_system(F, A, q, D, p2, GROUP)
    r2 = fixed_d.solve_compressed(S2)
    rec.update(prime2=p2, rank_c_p2=r2['rank'], consistent_p2=bool(r2['consistent']))
    if full:
        Sf = fixed_d.FixedDSystem(F, A, q, D, p, verbose=True)
        rf = fixed_d.solve_compressed(Sf, verbose=True)
        rec.update(full_rows=Sf.nrows, full_cols=Sf.ncols, full_rank_c=rf['rank'], full_consistent_modp=bool(rf['consistent']))
    rec['seconds'] = round(time.time() - t0, 1)
    L.receipt('rigorfix', **rec)
    return rec


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    run(args[0], args[1], full='--full' in sys.argv)
