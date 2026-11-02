"""Run the Koszul-reduced E2 system on a list of Jacobian denominators D = w(F) with rigor accounting.

Rigor: the system M x = b (x = (c, Q)) has an EXPLICIT rational structural kernel
   S = { (c(v), Q(v)) : v a vector field of degree d-5 } ,  N_v = 6 F v - v(F) E = sum_{i<j} c_ij(v) K_ij ,
   c_kj(v) = v_k Z_j - v_j Z_k (k<j),  Q(v) = 6 (D div v - v(D)),
   (plus, for d >= 9, second Koszul syzygies c_ij = sum_k e_ijk F_k which give N = 0).
If rank_p(M) == ncols - rank_p(S) then rank_Q(M) = rank_p(M) (sandwich: rank_p(M) <= rank_Q(M) <= ncols - dim_Q S
<= ncols - rank_p(S)), and modular inconsistency (rank_p([M|b]) = rank_p(M)+1 <= rank_Q([M|b])) is rigorous over Q.
"""
import sys, time, json, random
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl/vm')
from kred import *
from flint import nmod_mat
import numpy as np, scipy.sparse as sp

Z = [{tuple(1 if k == i else 0 for k in range(6)): 1} for i in range(6)]
EULER = [dict(z) for z in Z]

def struct_kernel_vectors(S, p):
    """explicit structural kernel vectors (mod p) in the column coordinates of KoszulSystem S."""
    d = S.d
    colidx = {c: j for j, c in enumerate(S.cols)}
    vecs = []
    D = S.D
    for k in range(6):
        for m in monomials(d - 5):
            v = [{} for _ in range(6)]; v[k] = {m: 1}
            # c_kj = v_k Z_j for j != k with sign: pair (min,max): c_{kj} K_{kj}, K_{kj} = F_j d_k - F_k d_j for k<j.
            # N_v = sum_j v_k Z_j K_{kj}  (for j<k, K_{kj} = -K_{jk} so c_{jk} = -v_k Z_j)
            vec = {}
            ok = True
            for j in range(6):
                if j == k: continue
                cpoly = mul(v[k], Z[j])
                for mm, cc in cpoly.items():
                    key = ('c', min(k, j), max(k, j), mm)
                    if key not in colidx: ok = False; break
                    sgn = 1 if k < j else -1
                    vec[colidx[key]] = (vec.get(colidx[key], 0) + sgn * cc) % p
                if not ok: break
            if not ok: continue   # parity-filtered system: this v not allowed
            Qv = scale(BD(D, v), 6)
            for mm, cc in Qv.items():
                key = ('Q', mm)
                if key not in colidx: ok = False; break
                vec[colidx[key]] = (vec.get(colidx[key], 0) + cc) % p
            if not ok: continue
            vecs.append(vec)
    if d >= 9:
        # second syzygies: c_ij = sum_k e_ijk F_k, e in S_{d-9} alternating -> N = 0, Q = 0
        for (a, b, c3) in [(a, b, c3) for a in range(6) for b in range(a + 1, 6) for c3 in range(b + 1, 6)]:
            for m in monomials(d - 9):
                vec = {}; ok = True
                # (d_2 e)_{ij} = sum_k e_{ijk} F_k with e_{abc}=m: c_ab += m F_c, c_ac -= m F_b, c_bc += m F_a
                for (i, j, k, sgn) in [(a, b, c3, 1), (a, c3, b, -1), (b, c3, a, 1)]:
                    cpoly = mul({m: 1}, reduce_mod(Fi_int[k], p))
                    for mm, cc in cpoly.items():
                        key = ('c', i, j, mm)
                        if key not in colidx: ok = False; break
                        vec[colidx[key]] = (vec.get(colidx[key], 0) + sgn * cc) % p
                    if not ok: break
                if ok: vecs.append(vec)
    return vecs

def struct_rank(S, p):
    vecs = struct_kernel_vectors(S, p)
    if not vecs: return 0, 0
    C = S.ncols
    flat = [0] * (len(vecs) * C)
    for r, vec in enumerate(vecs):
        for j, v in vec.items(): flat[r * C + j] = v
    M = nmod_mat(len(vecs), C, flat, p)
    return M.rank(), len(vecs)

def verify_struct_kernel(S, p, nsample=3):
    """check that sampled structural vectors are really in ker M (sanity of the formulas)."""
    vecs = struct_kernel_vectors(S, p)
    random.seed(0)
    for vec in random.sample(vecs, min(nsample, len(vecs))):
        acc = {}
        for j, v in vec.items():
            acc = add(acc, scale(S.colpolys[j], v, p), p)
        acc = {e: c for e, c in acc.items() if e in S.rows}
        if acc: return False
    return True

def run_one(name, w, p, parity_rows=None, compressed=False, verbose=True, family='', note=''):
    t0 = time.time()
    D = field_apply(w, Fint)
    if not D:
        return {'name': name, 'error': 'D = 0'}
    d = core.degree(D)
    pD = parity_of(D)
    pr = parity_rows if (parity_rows is not None and pD is not None) else None
    # coprimality: w(F) = hF iff w - hE/6 syzygy -> w in T; check F | D via divisibility test mod p
    divisible, _ = jacobian.divisible_by(Fint, D, p)
    S = KoszulSystem(w, parity_rows=pr, verbose=verbose)
    rec = {'exp': 'kred-family', 'family': family, 'name': name, 'd': d, 'parity_D': pD, 'rows_parity': pr,
           'F_divides_D_mod_p': bool(divisible), 'nrows': S.nrows, 'ncols': S.ncols, 'p': p, 'note': note, 'agent': 'excl-vm'}
    if compressed:
        M = sp.csr_matrix
        ri, ci, vals = [], [], []
        for j, col in enumerate(S.colpolys):
            for e, v in col.items():
                if e in S.rows: ri.append(S.rows[e]); ci.append(j); vals.append(v % p)
        for e, v in S.rhs.items():
            if e in S.rows: ri.append(S.rows[e]); ci.append(S.ncols); vals.append(v % p)
        Msp = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(S.nrows, S.ncols + 1))
        k = min(S.nrows, S.ncols + 1 + 40)
        Cm = fixed_d._compress(Msp, p, k, 1)
        Mn = nmod_mat(k, S.ncols + 1, Cm.ravel().tolist(), p)
        Mr, rank_aug = Mn.rref()
        piv = []
        for i in range(rank_aug):
            j = 0
            while j <= S.ncols and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        consistent = not (piv and piv[-1] == S.ncols)
        r = {'consistent': consistent, 'rank_p': rank_aug - (0 if consistent else 1), 'rank_aug_p': rank_aug}
        rec['solver'] = f'compressed k={k}'
    else:
        r = S.solve_mod(p)
        rec['solver'] = 'dense rref'
    rec.update({'consistent_p': r['consistent'], 'rank_p': r['rank_p'], 'rank_aug_p': r['rank_aug_p'], 'nullity_p': S.ncols - r['rank_p']})
    sr, nvec = struct_rank(S, p)
    rec['struct_kernel_rank_p'] = sr; rec['struct_kernel_vectors'] = nvec
    rec['struct_formulas_verified'] = verify_struct_kernel(S, p)
    rec['rank_matches_structural'] = bool(S.ncols - r['rank_p'] == sr)
    rec['rigorous_inconsistent_over_Q'] = bool((not r['consistent']) and rec['rank_matches_structural'] and not divisible)
    rec['seconds'] = round(time.time() - t0, 1)
    if r['consistent']:
        rec['note'] += ' !!! CONSISTENT mod p -> lift needed'
    if verbose: print(json.dumps(rec), flush=True)
    log(rec)
    return rec
