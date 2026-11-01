"""Reduced solver for denominators D in J(F): D = sum_i G_i F_i with G a (parity-compatible) polynomial vector field.
E1 is solved analytically: N = A*G + K(c), K(c) = sum_{i<j} c_ij (F_j e_i - F_i e_j).  E2: D div N - N(D) = F Q.
Unknowns: c_ij (degree deg D + 6q - 10), Q (degree 2 deg D - 6, Y-odd).  Returns consistency/rank mod p and exact lift."""
import sys, time, json, itertools, random
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import ctr
from ctr import *
from flint import nmod_mat
from fractions import Fraction
p0 = core.PRIMES[0]
Fi_exact = [core.pderiv(F, i) for i in range(6)]
Zlin = [P(v) for v in V]

def basis18():
    """Basis of J_6^even: X_j F_{X_i} (i,j<3) and Y_j F_{Y_i} (i,j>=3). Returns list of (name, G as list of 6 polys)."""
    out = []
    for i in range(3):
        for j in range(3):
            G = [{} for _ in range(6)]; G[i] = Zlin[j]
            out.append((f'{V[j]}*F_{V[i]}', G))
    for i in range(3, 6):
        for j in range(3, 6):
            G = [{} for _ in range(6)]; G[i] = Zlin[j]
            out.append((f'{V[j]}*F_{V[i]}', G))
    return out
B18 = basis18()

def D_of_G(G, p=None):
    D = {}
    for i in range(6):
        if G[i]:
            D = core.add(D, core.mul(G[i], Fi_exact[i] if p is None else core.reduce_mod(Fi_exact[i], p), p), p)
    return D

def koszul_system(G, p=p0, q=1):
    """Build the reduced E2 system mod p for D = G(F). Returns dict with consistent, rank, nullity, rows, cols, and solver data."""
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    Gp = [core.reduce_mod(g, p) for g in G]
    Dp = {}
    for i in range(6):
        if Gp[i]: Dp = core.add(Dp, core.mul(Gp[i], Fi[i], p), p)
    d = core.degree(Dp)
    if d < 0: return None
    Di = [core.pderiv(Dp, i, p) for i in range(6)]
    pars = fixed_d.numerator_parities(Dp, Ap)
    degc = d + 6*q - 10
    degQ = 2*d - 6
    # N0 = A*G
    N0 = [core.mul(Ap, Gp[i], p) if Gp[i] else {} for i in range(6)]
    def B_of(N):
        divN = {}; ND = {}
        for i in range(6):
            if N[i]:
                divN = core.add(divN, core.pderiv(N[i], i, p), p)
                ND = core.add(ND, core.mul(N[i], Di[i], p), p)
        return core.sub(core.mul(Dp, divN, p), ND, p)
    R0 = B_of(N0)
    cols = []; colpolys = []
    if degc >= 0:
        for i in range(6):
            for j in range(i+1, 6):
                parc = (pars[i] - core.parity_of(Fi[j]) if Fi[j] else 0) % 2 if Fi[j] else 0
                # parity of c: N_i component c*F_j must have parity pars[i]
                pj = 0 if j < 3 else 1
                parc = (pars[i] - pj) % 2
                for m in core.monomials(degc, parc):
                    cm = {m: 1}
                    N = [{} for _ in range(6)]
                    N[i] = core.mul(cm, Fi[j], p); N[j] = core.scale(core.mul(cm, Fi[i], p), p-1, p)
                    cols.append(('c', i, j, m)); colpolys.append(B_of(N))
    Fq = core.power(Fp, q, p)
    for m in core.monomials(degQ, core.parity_of(A)) if degQ >= 0 else []:
        cols.append(('Q', -1, -1, m)); colpolys.append(core.scale(core.shift(Fq, m, 1, p), p-1, p))
    rhs = core.scale(R0, p-1, p)
    rows = {}
    for c in colpolys:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * (C + 1))
    for j, c in enumerate(colpolys):
        for e, v in c.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(R, C + 1, flat, p)
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    # rank of M without rhs
    M0 = nmod_mat(R, C, [flat[r*(C+1)+j] for r in range(R) for j in range(C)], p)
    rank0 = M0.rank()
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
    return {'consistent': consistent, 'rank': rank0, 'rank_aug': rank, 'nullity': C - rank0, 'rows': R, 'cols': C,
            'd': d, 'seconds': round(time.time() - t0, 2), 'x': x, 'cols_meta': cols, 'p': p}

def lift_solution(G, primes=None, q=1, verbose=False):
    """CRT + rational reconstruction over several primes; exact verification."""
    primes = primes or core.PRIMES
    acc = None; meta = None
    for k, p in enumerate(primes[:8]):
        r = koszul_system(G, p=p, q=q)
        if not r['consistent']: return None
        x = r['x']
        if meta is None: meta = r['cols_meta']
        if acc is None: acc = (x, p)
        else:
            xs, m = acc
            xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, x)]
            acc = (xs, m*p)
        xs, m = acc
        rat = [core.ratrecon(a, m) for a in xs]
        if all(v is not None for v in rat):
            # assemble N, Q exactly
            D = D_of_G(G)
            N = [core.mul(A, G[i]) if G[i] else {} for i in range(6)]
            Q = {}
            for (kind, i, j, mon), v in zip(meta, rat):
                if v == 0: continue
                if kind == 'c':
                    cm = {mon: v}
                    N[i] = core.add(N[i], core.mul(cm, Fi_exact[j]))
                    N[j] = core.sub(N[j], core.mul(cm, Fi_exact[i]))
                else:
                    Q[mon] = v
            ok, info = core.verify_witness(F, A, q, D, N, Q)
            if ok: return {'D': D, 'N': N, 'Q': Q, 'primes': k+1}
            elif verbose: print('lift verify failed, continue', file=sys.stderr)
    return None

def test_G(G, name, family='J6', note=''):
    r = koszul_system(G)
    D = D_of_G(G)
    ent = {'engine': 2, 'D': name[:400], 'd': r['d'], 'family': family, 'e1': True, 'consistent': bool(r['consistent']),
           'rank': r['rank'], 'nullity': r['nullity'], 'seconds': r['seconds'],
           'note': (note + f" | reduced Koszul system rows={r['rows']} cols={r['cols']} rank_aug={r['rank_aug']} (D in J(F), E1 automatic)").strip(' |')}
    if r['consistent']:
        sol = lift_solution(G)
        if sol:
            ent['solved'] = True; ent['verified'] = True
            ent['witness'] = {'d': r['d'], 'D': core.to_str(sol['D']), 'N': [core.to_str(n) for n in sol['N']], 'Q': core.to_str(sol['Q'])}
            with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
                json.dump({'instanceId': IID, 'family': family, 'name': name, **ent['witness']}, f, indent=1)
            print('\n*** SOLVED (J-family) ***', name, '\n', flush=True)
        else:
            ent['note'] += ' | CONSISTENT mod p but exact lift failed'
    log(ent)
    print(f"d={ent['d']} cons={str(ent['consistent']):5s} rank={ent['rank']:4d} null={ent['nullity']:4d} rank_aug={r['rank_aug']:4d} {ent['seconds']:5.1f}s {family:10s} {name[:80]}", flush=True)
    return ent, r

def G_from_coeffs(s):
    G = [{} for _ in range(6)]
    for sk, (name, Gk) in zip(s, B18):
        if sk:
            for i in range(6):
                if Gk[i]: G[i] = core.add(G[i], core.scale(Gk[i], Fraction(sk)))
    return G

if __name__ == '__main__':
    # cross-check against the general solver on one element
    name, G = B18[0]
    ent, r = test_G(G, name, note='crosscheck')
    D = D_of_G(G)
    out = fixed_d.solve_fixed_D_fast(F, A, 1, D, want_exact=False)
    print('general solver on same D:', {k: out[k] for k in ('consistent', 'rank', 'nullity', 'rows', 'cols')}, flush=True)
