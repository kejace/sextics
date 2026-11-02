"""Served target spaces (sum principle). For a G-semi-invariant denominator P (q = 1), compute
   served(P)^{chi} = { [B] in (S_6^{chi} / J_6^{chi}) : B Omega/F^2 has a certificate with denominator P }   (mod p),
for the sector chi = chi_A of the target. At q = 1 served(P1 P2) = served(P1) + served(P2) mod J, so [A] in the span of
several served spaces gives an explicit product denominator. usage: python3 w2_served.py cands.jsonl [p] [procs]
Each candidate line: {"D": expr, "name": ..., "family": ...}. Writes w2-served-<name-hash>.json and logs one line per candidate
plus a final span test."""
import os, sys, time, json, hashlib
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import multiprocessing as mp
from w2core import *
from flint import nmod_mat

F, A, q = load()
G = full_group(F, A)
S = Sect(F, A, G)
p = int(sys.argv[2]) if len(sys.argv) > 2 else PRIMES[0]
procs = int(sys.argv[3]) if len(sys.argv) > 3 else 1
Fz = sect.integral(F); dFz = [sect.ideriv(Fz, i) for i in range(6)]
chiA = S.chiA
Tm = S.mons(6, chiA)                      # 41 target monomials in the sector of A
Tidx = {m: k for k, m in enumerate(Tm)}
# J_6^{chiA} (mod p) reducer rows for the quotient
Jgens6 = [sect.imul({m: 1}, dFz[i]) for i in range(6) for m in S.mons(1, S.chiJgen(chiA, i))]
def to_vec(P):
    v = [0] * len(Tm)
    for e, c in P.items():
        v[Tidx[e]] = c % p
    return v
Jmat = nmod_mat(len(Jgens6), len(Tm), [x for g in Jgens6 for x in to_vec(g)], p) if Jgens6 else None
rankJ = Jmat.rank() if Jmat is not None else 0
Avec = to_vec(sect.integral(A))

def served(D):
    """Basis (mod p) of {B in S_6^{chiA}: B D in J, with full E2} i.e. targets served by D, as vectors over Tm; plus mod-J rank."""
    Dz = sect.integral(D)
    c = sect.chi_of(Dz, S.chi)
    d = sum(next(iter(Dz)))
    chiN = [S.chiN(c, i) for i in range(6)]
    # columns: N_i monomials (deg d+1), Q monomials (deg 2d-6, chi = chiA + 2c), target monomials B (sector chiA): equations
    # E1: sum N_i Fz_i - B Dz = 0 (rows deg d+6, chi chiA + c); E2: Dz div N - N(Dz) - Fz Q = 0 (rows deg 2d, chi chiA + 2c)
    r1 = {m: k for k, m in enumerate(S.mons(d + 6, sect.cadd(chiA, c)))}
    r2 = {m: k for k, m in enumerate(S.mons(2 * d, S.chiQ(c)))}
    n1 = len(r1); R = n1 + len(r2)
    cols = []; meta = []
    dDz = [sect.ideriv(Dz, i) for i in range(6)]
    for i in range(6):
        for m in S.mons(d + 1, chiN[i]):
            col = {}
            for e, v in sect.imul({m: 1}, dFz[i]).items():
                col[r1[e]] = (col.get(r1[e], 0) + v) % p
            t2 = sect.imul(Dz, sect.ideriv({m: 1}, i))
            for e, v in sect.imul({m: 1}, dDz[i]).items():
                t2[e] = t2.get(e, 0) - v
            for e, v in t2.items():
                if v % p:
                    k = n1 + r2[e]; col[k] = (col.get(k, 0) + v) % p
            cols.append(col); meta.append('N')
    for m in (S.mons(2 * d - 6, S.chiQ(c)) if 2 * d - 6 >= 0 else []):
        col = {}
        for e, v in sect.imul({m: 1}, Fz).items():
            k = n1 + r2[e]; col[k] = (col.get(k, 0) - v) % p
        cols.append(col); meta.append('Q')
    nNQ = len(cols)
    for m in Tm:
        col = {}
        for e, v in sect.imul({m: 1}, Dz).items():
            col[r1[e]] = (col.get(r1[e], 0) - v) % p
        cols.append(col); meta.append('B')
    C = len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(cols):
        for r_, v in col.items():
            flat[r_ * C + j] = v
    M = nmod_mat(R, C, flat, p)
    X, nullity = M.nullspace()
    # projection of the nullspace onto the B coordinates
    if nullity == 0:
        return [], 0, {'rows': R, 'cols': C, 'nullity': 0}
    proj = nmod_mat(nullity, len(Tm), [int(X[nNQ + k, j]) for j in range(nullity) for k in range(len(Tm))], p)
    Pr, rk = proj.rref()
    basis = [[int(Pr[i, k]) for k in range(len(Tm))] for i in range(rk)]
    # rank modulo J_6
    if Jmat is not None and rk:
        stack = nmod_mat(len(Jgens6) + rk, len(Tm), [x for g in Jgens6 for x in to_vec(g)] + [x for b in basis for x in b], p)
        rank_modJ = stack.rank() - rankJ
    else:
        rank_modJ = rk
    return basis, rank_modJ, {'rows': R, 'cols': C, 'nullity': nullity, 'raw_rank': rk}

def work(cand):
    t0 = time.time()
    D = core.parse(cand['D'])
    rec = {'batch': 'served', 'name': cand.get('name', cand['D'][:60]), 'family': cand.get('family'), 'd': core.degree(D), 'p': p}
    try:
        Dz = sect.integral(D)
        c = sect.chi_of(Dz, S.chi)
        rec['sector'] = list(c)
        basis, rank_modJ, info = served(D)
        rec.update(info); rec['served_dim_modJ'] = rank_modJ; rec['served_dim_raw'] = len(basis)
        # does the served space contain [A] (mod J)?
        if basis:
            stack = nmod_mat(len(Jgens6) + len(basis) + 1, len(Tm), [x for g in Jgens6 for x in to_vec(g)] + [x for b in basis for x in b] + Avec, p)
            rec['contains_A'] = (stack.rank() == rankJ + rank_modJ)
        else:
            rec['contains_A'] = False
        rec['basis'] = basis
    except Exception as ex:
        import traceback
        rec['error'] = repr(ex); rec['tb'] = traceback.format_exc()[-500:]
    rec['seconds'] = round(time.time() - t0, 1)
    return rec

if __name__ == '__main__':
    cands = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
    print(f'{len(cands)} candidates; target sector dim {len(Tm)}, J_6 part rank {rankJ}, R_6 sector dim {len(Tm) - rankJ}', flush=True)
    allb = []
    with mp.get_context('fork').Pool(procs) as pool:
        for rec in pool.imap_unordered(work, cands):
            basis = rec.pop('basis', [])
            if basis:
                allb.extend(basis)
                with open(f'{WD}/w2-served-{hashlib.md5(rec["name"].encode()).hexdigest()[:10]}.json', 'w') as f:
                    json.dump({'name': rec['name'], 'basis': basis, 'Tm': [core.to_str({m: 1}) for m in Tm]}, f)
            w2log(rec)
            print(f"d={rec.get('d')} served_modJ={rec.get('served_dim_modJ')} containsA={rec.get('contains_A')} {rec.get('seconds')}s {rec.get('family')} {rec['name'][:60]} {rec.get('error', '')}", flush=True)
    # span test
    if allb:
        stack = nmod_mat(len(Jgens6) + len(allb), len(Tm), [x for g in Jgens6 for x in to_vec(g)] + [x for b in allb for x in b], p)
        rk = stack.rank()
        stack2 = nmod_mat(len(Jgens6) + len(allb) + 1, len(Tm), [x for g in Jgens6 for x in to_vec(g)] + [x for b in allb for x in b] + Avec, p)
        inA = stack2.rank() == rk
        res = {'batch': 'served-span', 'n_candidates': len(cands), 'span_dim_modJ': rk - rankJ, 'sector_dim_R6': len(Tm) - rankJ, 'A_in_span': inA, 'p': p}
    else:
        res = {'batch': 'served-span', 'n_candidates': len(cands), 'span_dim_modJ': 0, 'sector_dim_R6': len(Tm) - rankJ, 'A_in_span': False, 'p': p}
    w2log(res); print(json.dumps(res), flush=True)
