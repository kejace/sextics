"""Koszul correction-map fixed-D solver for klein-coupled-character (q = 1).

Given D and an E1 lift N0 with N0(F) = A D, every numerator is N = N0 + K(p), K(p)_i = sum_j p_ij F_j
(p antisymmetric; F_0..F_5 regular sequence => all syzygies Koszul).  Then E = F * B_D(N) and the
certificate condition is B_D(N0) + sum_{i<j} p_ij L_ij = F Q with
  B_D(N) = D div N - N(D),   L_ij(m) = D[(d_i m)F_j - (d_j m)F_i] - m(F_j D_i - F_i D_j).
Unknowns: p_ij monomials (deg d-4, parity pA+pD+pi_i+pi_j, char a+k+w_i+w_j) and Q (deg 2d-6, parity pA, char a+2k).
"""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, yparity
from kcc import W, wt, sigma_weight, wmonomials, log, P0
from flint import nmod_mat
from fractions import Fraction

PI = (0, 0, 0, 1, 1, 1)
PAIRS = [(i, j) for i in range(NV) for j in range(i+1, NV)]

def dense_from_cols(colpolys, rhs, p, rows=None):
    if rows is None:
        rows = {}
        for col in colpolys:
            for e in col: rows.setdefault(e, len(rows))
        if rhs:
            for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(colpolys)
    Wd = C + (1 if rhs is not None else 0)
    flat = [0] * (R * Wd)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rows[e] * Wd + j] = v % p
    if rhs is not None:
        for e, v in rhs.items():
            flat[rows[e] * Wd + C] = v % p
    return nmod_mat(R, Wd, flat, p), rows

def rref_solve(M, C, p):
    """M = [A | b] (R x (C+1)). Returns (consistent, rank, pivots, x with free vars = 0)."""
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
    return consistent, rank, pivots, x

class Ctx:
    def __init__(self, F, A, p):
        self.p = p
        self.F = core.reduce_mod(F, p); self.A = core.reduce_mod(A, p)
        self.Fi = [pderiv(self.F, i, p) for i in range(NV)]
        self.pA = core.parity_of(A); self.a = sigma_weight(A)

def e1_lift(ctx, D, weights=True):
    """Solve N(F) = A D (mod p), N_i of degree d+1. Returns (ok, N, pivots)."""
    p = ctx.p; D = core.reduce_mod(D, p)
    d = core.degree(D); pD = core.parity_of(D); k = sigma_weight(D)
    restrict = weights and k is not None
    cols, colpolys = [], []
    for i in range(NV):
        par = (ctx.pA - pD - PI[i]) % 2
        kN = (ctx.a + k + W[i]) % 7 if restrict else None
        for m in wmonomials(d + 1, par, kN):
            cols.append((i, m)); colpolys.append(shift(ctx.Fi[i], m, 1, p))
    rhs = mul(ctx.A, D, p)
    M, rows = dense_from_cols(colpolys, rhs, p)
    ok, rank, pivots, x = rref_solve(M, len(cols), p)
    if not ok: return False, None, pivots
    N = [{} for _ in range(NV)]
    for (i, m), v in zip(cols, x):
        if v: N[i][m] = v
    return True, N, pivots

def koszul_columns(ctx, D, N0, weights=True):
    p = ctx.p; D = core.reduce_mod(D, p)
    d = core.degree(D); pD = core.parity_of(D); k = sigma_weight(D)
    restrict = weights and k is not None
    Di = [pderiv(D, i, p) for i in range(NV)]
    Fi = ctx.Fi; F = ctx.F
    def B(N):  # D div N - N(D)
        acc = {}
        for i in range(NV):
            if N[i]:
                acc = add(acc, mul(D, pderiv(N[i], i, p), p), p)
                acc = sub(acc, mul(N[i], Di[i], p), p)
        return acc
    cols, colpolys = [], []
    if d >= 4:
        for (i, j) in PAIRS:
            par = (ctx.pA + pD + PI[i] + PI[j]) % 2
            kp = (ctx.a + k + W[i] + W[j]) % 7 if restrict else None
            # L_ij(m) = D[(d_i m)F_j - (d_j m)F_i] - m(F_j D_i - F_i D_j)
            G = sub(mul(Fi[j], Di[i], p), mul(Fi[i], Di[j], p), p)
            DFj = mul(D, Fi[j], p); DFi = mul(D, Fi[i], p)
            for m in wmonomials(d - 4, par, kp):
                col = scale(shift(G, m, 1, p), p - 1, p)
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(DFj, tuple(m2), m[i], p), p)
                if m[j]:
                    m2 = list(m); m2[j] -= 1
                    col = sub(col, shift(DFi, tuple(m2), m[j], p), p)
                cols.append(('p', i, j, m)); colpolys.append(col)
    nP = len(cols)
    if 2 * d - 6 >= 0:
        kQ = (ctx.a + 2 * k) % 7 if restrict else None
        for m in wmonomials(2 * d - 6, ctx.pA, kQ):
            cols.append(('Q', m)); colpolys.append(scale(shift(F, m, 1, p), p - 1, p))
    rhs = scale(B(N0), p - 1, p)   # move B(N0) to rhs:  sum p L + (-F Q) = -B(N0)
    return cols, colpolys, rhs, nP

def koszul_test(ctx, D, N0, weights=True, verbose=False):
    t0 = time.time()
    cols, colpolys, rhs, nP = koszul_columns(ctx, D, N0, weights)
    M, rows = dense_from_cols(colpolys, rhs, ctx.p)
    ok, rank, pivots, x = rref_solve(M, len(cols), ctx.p)
    rank0 = rank if ok else rank - 1   # rank of the coefficient matrix (augmented rank minus 1 when inconsistent)
    res = {'consistent': ok, 'rank': rank0, 'nullity': len(cols) - rank0, 'rows': len(rows), 'cols': len(cols),
           'nP': nP, 'nQ': len(cols) - nP, 'seconds': round(time.time() - t0, 2), 'pivots': pivots, 'x': x, 'colkeys': cols}
    if verbose:
        print(f"[koszul p={ctx.p}] d={core.degree(D)} k={sigma_weight(D)} rows={res['rows']} cols={res['cols']} (p {nP}, Q {res['nQ']}) "
              f"consistent={ok} rank={rank} nullity={res['nullity']} {res['seconds']}s", file=sys.stderr)
    return res

def assemble(ctx, N0, res):
    """N = N0 + K(p), Q from a solution vector (mod p)."""
    p = ctx.p
    N = [dict(n) for n in N0]; Q = {}
    for key, v in zip(res['colkeys'], res['x']):
        if not v: continue
        if key[0] == 'p':
            _, i, j, m = key
            N[i] = add(N[i], shift(ctx.Fi[j], m, v, p), p)
            N[j] = sub(N[j], shift(ctx.Fi[i], m, v, p), p)
        else:
            Q[key[1]] = v
    return N, Q

def full_test(F, A, D, name, family, note='', weights=True, verbose=False, p=P0, want_log=True):
    """E1 lift + Koszul test for a rational (or mod-p) D. Logs. Returns entry."""
    t0 = time.time()
    ctx = Ctx(F, A, p)
    d = core.degree(D); k = sigma_weight(D)
    ok, N0, piv = e1_lift(ctx, D, weights)
    entry = {'D': name if len(name) < 300 else name[:300] + '...', 'd': d, 'family': family, 'k': k, 'e1': ok, 'engine': 'koszul',
             'consistent': None, 'rank': None, 'nullity': None, 'note': note}
    if ok:
        res = koszul_test(ctx, D, N0, weights, verbose)
        entry.update({'consistent': res['consistent'], 'rank': res['rank'], 'nullity': res['nullity'], 'rows': res['rows'], 'cols': res['cols']})
        entry['_res'] = res; entry['_N0'] = N0; entry['_ctx'] = ctx
    entry['seconds'] = round(time.time() - t0, 2)
    if want_log: log({kk: v for kk, v in entry.items() if not kk.startswith('_')})
    if verbose and not ok: print(f'[e1] FAIL d={d} k={k} {name[:60]}', file=sys.stderr)
    return entry

# ---------------------------------------------------------------- rational lifting
def lift_vector_multiprime(solver, nprimes=8, verbose=False):
    """solver(p) -> (pivots, x list of ints mod p) canonical given pivots; CRT + ratrecon."""
    acc = None; piv0 = None
    for kk, p in enumerate(core.PRIMES[:nprimes]):
        piv, x = solver(p)
        if x is None: return None
        if piv0 is None or piv != piv0:
            if piv0 is not None and verbose: print('[lift] pivot mismatch, restarting', file=sys.stderr)
            piv0 = piv; acc = (list(x), p); continue
        xs, m = acc
        xs = [core.crt_pair(aa, m, b, p)[0] for aa, b in zip(xs, x)]
        acc = (xs, m * p)
        rat = [core.ratrecon(aa, acc[1]) for aa in acc[0]]
        if all(v is not None for v in rat):
            return rat
    return None

def lift_certificate(F, A, D, weights=True, nprimes=8, verbose=False):
    """D rational. Lift N0 and the Koszul solution to Q, verify with core.verify_witness."""
    def solver(p):
        ctx = Ctx(F, A, p)
        ok, N0, piv1 = e1_lift(ctx, D, weights)
        if not ok: return None, None
        res = koszul_test(ctx, D, N0, weights)
        if not res['consistent']: return None, None
        N, Q = assemble(ctx, N0, res)
        # serialize in a canonical order: all N monomials then Q monomials over fixed monomial lists
        d = core.degree(D); pD = core.parity_of(D)
        vec = []
        for i in range(NV):
            for m in monomials(d + 1, (ctx.pA - pD - PI[i]) % 2): vec.append(N[i].get(m, 0))
        for m in monomials(2 * d - 6, ctx.pA): vec.append(Q.get(m, 0))
        return (tuple(piv1), tuple(res['pivots'])), vec
    rat = lift_vector_multiprime(solver, nprimes, verbose)
    if rat is None: return None
    d = core.degree(D); pD = core.parity_of(D); pA = core.parity_of(A)
    N = [{} for _ in range(NV)]; Q = {}; idx = 0
    for i in range(NV):
        for m in monomials(d + 1, (pA - pD - PI[i]) % 2):
            if rat[idx]: N[i][m] = rat[idx]
            idx += 1
    for m in monomials(2 * d - 6, pA):
        if rat[idx]: Q[m] = rat[idx]
        idx += 1
    ok, info = core.verify_witness(F, A, 1, D, N, Q)
    return {'ok': ok, 'info': info, 'N': N, 'Q': Q}

# ---------------------------------------------------------------- E1 kernel (canonical basis) and J blocks
def jac_block_cols(ctx, deg, par, k):
    """columns F_i * m spanning J_deg in the (par, k) block; returns (keys, polys)."""
    p = ctx.p; keys, polys = [], []
    for i in range(NV):
        pm = (par + PI[i]) % 2 if par is not None else None
        km = (k + W[i]) % 7 if k is not None else None
        for m in wmonomials(deg - 5, pm, km):
            keys.append((i, m)); polys.append(shift(ctx.Fi[i], m, 1, p))
    return keys, polys

def e1_kernel(ctx, d, k, parD=0, verbose=False):
    """Canonical basis (mod p) of {D in S_d^(parD,k) : A D in J} modulo J_d, as a list of polys,
    plus a basis of J_d block (keys, polys). Kernel reps are reduced against J_d pivots (canonical)."""
    p = ctx.p
    Dm = wmonomials(d, parD, k)
    if not Dm: return [], ([], [])
    tgt = wmonomials(d + 6, (parD + ctx.pA) % 2, (k + ctx.a) % 7 if k is not None else None)
    rows = {m: i for i, m in enumerate(tgt)}
    jkeys, jpolys = jac_block_cols(ctx, d + 6, (parD + ctx.pA) % 2, (k + ctx.a) % 7 if k is not None else None)
    # left null space L of J block: functionals on S_{d+6} vanishing on J
    MJ, _ = dense_from_cols(jpolys, None, p, rows) if jpolys else (nmod_mat(len(rows), 0, p), None)
    L, nl = MJ.transpose().nullspace()   # columns of L (R x nl) are functionals (only first nl columns meaningful)
    # multiplication matrix: for each D monomial, A*m reduced by functionals -> nl-vector
    AD, _ = dense_from_cols([shift(ctx.A, m, 1, p) for m in Dm], None, p, rows)   # R x nD
    Lt = nmod_mat(nl, len(rows), [int(L[r, c]) for c in range(nl) for r in range(len(rows))], p)
    Mult = Lt * AD   # nl x nD
    K, nk = Mult.nullspace()   # kernel in D-coordinates (canonical from rref)
    kervecs = [[int(K[r, c]) for r in range(len(Dm))] for c in range(nk)]
    # J_d block in D-space
    jdkeys, jdpolys = jac_block_cols(ctx, d, parD, k)
    jd = 0; Jmat = None
    if jdpolys:
        Jmat, _ = dense_from_cols(jdpolys, None, p, {m: i for i, m in enumerate(Dm)})
        jd = Jmat.rank()
    # reduce kernel modulo J_d: stack [J_d cols | kernel vecs], rref, keep kernel vecs independent mod J
    reps = []
    if nk:
        nD = len(Dm)
        flat = []
        allcols = ([[int(Jmat[r, c]) for r in range(nD)] for c in range(len(jdpolys))] if jdpolys else []) + kervecs
        C = len(allcols)
        flat = [0] * (nD * C)
        for j, col in enumerate(allcols):
            for r, v in enumerate(col):
                if v: flat[r * C + j] = v
        Mall = nmod_mat(nD, C, flat, p)
        Mr, rk = Mall.rref()
        piv = []
        for i in range(rk):
            j = 0
            while j < C and int(Mr[i, j]) == 0: j += 1
            piv.append(j)
        nJ = len(jdpolys)
        chosen = [j - nJ for j in piv if j >= nJ]
        for c in chosen:
            reps.append({Dm[r]: v for r, v in enumerate(kervecs[c]) if v})
    if verbose:
        print(f'[e1-kernel] d={d} k={k} parD={parD}: #Dmon={len(Dm)} dimR_(d+6)={nl} kernel(total)={nk} J_d={jd} nontrivial={len(reps)}', file=sys.stderr)
    return reps, (jdkeys, jdpolys)

def load():
    F, A, q, rec = core.load_instance('klein-coupled-character')
    return F, A, q
