"""Normalized E2 test (mod p) for an E1-admissible denominator D with an exact E1 numerator N0 (N0(F) = A D):
   exists K skew (15 entries, degree d-4) and C in S_{2d-6} with  B_D(N0) + B_D(K grad F) = F C ?
   (lossless: the Euler gauge g*Z drops out of E2 identically; Koszul syzygies exhaust the E1 solution set
    because grad F is a regular sequence).  Inconsistent mod p => inconsistent over Q when the mod-p rank of the
    homogeneous block equals its Q-rank; we report the modular statement and (optionally) a second prime.
Usage: python3 e2norm.py pairs.json [indices...]   |   python3 e2norm.py --selftest
       python3 e2norm.py --D "<expr>" [--F "<expr>"]  (D in J: N0 := A w is built automatically if D = sum w_i F_i given as --w "e0|e1|...|e5")
"""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/wave2-excl')
from qx import *

def parse_exp(dct): return {tuple(int(x) for x in k.split()): v for k, v in dct.items()}

def e2_normalized(F, A, D, N0, p, tag='', family='e2norm', do_log=True):
    t0 = time.time()
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p); D = core.reduce_mod(D, p) if any(isinstance(v, Fraction) for v in D.values()) else D
    N0 = [core.reduce_mod(n, p) if any(isinstance(v, Fraction) for v in n.values()) else n for n in N0]
    Fi = [pderiv(Fp, i, p) for i in range(NV)]
    d = core.degree(D)
    lhs = {}
    for i in range(NV): lhs = add(lhs, mul(N0[i], Fi[i], p), p)
    assert not sub(lhs, mul(Ap, D, p), p), 'E1 fails for the given (D, N0)'
    Di = [pderiv(D, i, p) for i in range(NV)]
    divN = {}
    for i in range(NV): divN = add(divN, pderiv(N0[i], i, p), p)
    ND = {}
    for i in range(NV): ND = add(ND, mul(N0[i], Di[i], p), p)
    rhs = scale(sub(mul(D, divN, p), ND, p), p - 1, p)      # -B_D(N0)
    pars = fixed_d.numerator_parities(D, Ap)
    DFij = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            DFij[(i, j)] = sub(mul(Di[i], Fi[j], p), mul(Di[j], Fi[i], p), p)
    cols = []; names = []
    for i in range(NV):
        for j in range(i + 1, NV):
            pj = 0 if j < 3 else 1
            for m in monomials(d - 4, (pars[i] - pj) % 2):
                mm = {m: 1}
                # N_i += K_ij F_j, N_j -= K_ij F_i ;  B_D contribution = D*(d_i(K F_j) + d_j(-K F_i)) - (K F_j D_i - K F_i D_j)
                t = mul(D, sub(mul(pderiv(mm, i, p), Fi[j], p), mul(pderiv(mm, j, p), Fi[i], p), p), p)
                t = sub(t, mul(mm, DFij[(i, j)], p), p)
                cols.append(t); names.append(('K', i, j, m))
    parC = (core.parity_of(Ap) + 2 * core.parity_of(D)) % 2 if False else None
    # parity of C: B_D(N0) has parity of A*D^2 = parity(A) (= odd) when D pure; F even => C odd
    pC = (core.parity_of(Ap)) % 2
    for m in monomials(2 * d - 6, pC):
        cols.append(scale(shift(Fp, m, 1, p), p - 1, p)); names.append(('C', m))
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    Rn, C = len(rows), len(cols)
    flat = [0] * (Rn * (C + 1))
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items(): flat[rows[e] * (C + 1) + C] = v
    M = nmod_mat(Rn, C + 1, flat, p); del flat
    t1 = time.time()
    r_aug = M.rank()
    M0 = nmod_mat(Rn, C, [0] * (Rn * C), p)
    flat0 = [0] * (Rn * C)
    for j, c in enumerate(cols):
        for e, v in c.items(): flat0[rows[e] * C + j] = v
    M0 = nmod_mat(Rn, C, flat0, p); del flat0
    r0 = M0.rank()
    rec = {'exp': 'e2_normalized', 'tag': tag, 'd': d, 'p': p, 'family': family, 'e1': True,
           'consistent_mod_p': (r0 == r_aug), 'rank': r0, 'rank_aug': r_aug, 'rows': Rn, 'cols': C, 'nullity': C - r0,
           'terms_D': len(D), 'seconds': round(time.time() - t0, 1)}
    if do_log: log(rec)
    print(json.dumps(rec), flush=True)
    return rec

if __name__ == '__main__':
    F, A, q = load()
    if sys.argv[1] == '--selftest':
        # Fermat control: D = X0^6+X1^6+X2^6, N0 = A E_X / 6 (N0(F0) = A D) must be E2-consistent; and D' = X0^6+X1^6+Y0^6?
        F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
        D0 = core.parse('X0**6+X1**6+X2**6')
        N0 = [core.parse('X0**4*Y0**3/6'), core.parse('X0**3*X1*Y0**3/6'), core.parse('X0**3*X2*Y0**3/6'), {}, {}, {}]
        r = e2_normalized(F0, A, D0, N0, P2, tag='fermat-control-join', family='selftest', do_log=False)
        assert r['consistent_mod_p'], 'selftest 1 failed (must be consistent)'
        # control 2: on Fermat, D = Y0 F_0/... use D = X0*F0_0 (polar, in J), N0 = A*(X0,0,0,0,0,0): expected inconsistent? (unknown) just report
        D1 = core.parse('X0**6+X1**6+Y0**6')   # join of a different type on Fermat: E1 holds? A*D1 in J? check
        ok, info = fixed_d.e1_solvable(F0, A, 1, D1, p=P2)
        print('control 2: E1 for X0^6+X1^6+Y0^6 on Fermat:', ok)
        print('SELFTEST OK')
        sys.exit(0)
    fn = sys.argv[1]; which = [int(x) for x in sys.argv[2:] if x.isdigit()]
    data = json.load(open(fn)); p = data['p']; d = data['d']
    for ri, rep in enumerate(data['pairs']):
        if which and ri not in which: continue
        D = parse_exp(rep['D']); N0 = [parse_exp(x) for x in rep['N']]
        e2_normalized(F, A, D, N0, p, tag=f'{os.path.basename(fn)}#{ri}', family=f'e1kernel-d{d}-E2norm')
