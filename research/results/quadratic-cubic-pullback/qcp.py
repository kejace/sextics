"""Helpers for the quadratic-cubic-pullback attack (work dir module).

Usage:  import sys; sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback'); import qcp
"""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction
import itertools

IID = 'quadratic-cubic-pullback'
WD = '/tmp/hodge/results/quadratic-cubic-pullback'
LOG = os.path.join(WD, 'log.jsonl')

F, A, q, rec = core.load_instance(IID)

QSTR = [
    "X0**2 + 5*X1*X2/7 + 5*Y0*Y1/7",
    "5*X0*X2/11 + X1**2 + 5*Y1*Y2/11",
    "5*X0*X1/13 + X2**2 + 5*Y0*Y2/13",
    "5*X0*X1/17 + Y0**2 + 5*Y1*Y2/17",
    "5*X1*X2/19 + 5*Y0*Y2/19 + Y1**2",
    "5*X0*X2/23 + 5*Y0*Y1/23 + Y2**2",
]
Q = [core.parse(s) for s in QSTR]

FERMAT = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')

P = core.parse


def det(M, p=None):
    """Determinant of a matrix of polynomials (dict polys) by Laplace/Bareiss-free expansion (n<=6)."""
    n = len(M)
    if n == 1:
        return M[0][0]
    out = {}
    for j in range(n):
        if not M[0][j]:
            continue
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        term = core.mul(M[0][j], det(minor, p), p)
        if j % 2:
            term = core.scale(term, -1, p)
        out = core.add(out, term, p)
    return out


def jac_det(Phi, p=None):
    M = [[core.pderiv(Phi[i], j, p) for j in range(6)] for i in range(6)]
    return det(M, p)


def log(D, family, e1, consistent, rank, nullity, seconds, note='', d=None, extra=None):
    rec_ = {"D": core.to_str(D) if isinstance(D, dict) else str(D),
            "d": core.degree(D) if isinstance(D, dict) else d,
            "family": family, "e1": e1, "consistent": consistent,
            "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    if extra:
        rec_.update(extra)
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec_) + '\n')
    return rec_


def test(D, family, note='', F_=None, A_=None, q_=None, full=True, fast=True, verbose=False, exact=True):
    """E1 test then (optionally) full fixed-D solve.  Logs the result.  Returns dict."""
    F_ = F if F_ is None else F_
    A_ = A if A_ is None else A_
    q_ = q if q_ is None else q_
    t0 = time.time()
    if core.parity_of(D) is None:
        return log(D, family, None, None, -1, -1, 0, note + ' [mixed parity, skipped]')
    ok, info = fixed_d.e1_solvable(F_, A_, q_, D)
    res = {'e1': ok, 'e1info': info}
    if not ok or not full:
        r = log(D, family, ok, None if ok else False, info['rank'], info['cols'] - info['rank'],
                time.time() - t0, note + (' [E1 fails]' if not ok else ' [E1 only]'))
        res.update(r)
        return res
    if fast:
        out = fixed_d.solve_fixed_D_fast(F_, A_, q_, D, verbose=verbose, want_exact=exact)
        cons, rank, nullity = out['consistent'], out['rank'], out['nullity']
        res.update(out)
        if cons and exact and out.get('exact'):
            okv, vinfo = core.verify_witness(F_, A_, q_, D, out['N'], out['Q'])
            res['verified'] = okv
            note = note + f' [SOLVED exact verify={okv}]'
    else:
        status, data = fixed_d.solve_fixed_D(F_, A_, q_, D, verbose=verbose)
        cons = status == 'solved'
        rank, nullity = data.get('rank', -1), data.get('nullity', -1)
        res.update({'status': status, 'data': data})
        if cons:
            res['N'], res['Q'] = data['N'], data['Q']
            res['verified'] = True
            note = note + ' [SOLVED exact]'
    r = log(D, family, ok, cons, rank, nullity, time.time() - t0, note)
    res.update(r)
    return res


def image_dim(D, F_=None, q_=None, p=None):
    F_ = F if F_ is None else F_
    q_ = q if q_ is None else q_
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F_, q_, D, p=p)
    return dim, basis, S, time.time() - t0


def save_witness(D, N, Qp, path=os.path.join(WD, 'witness-candidate.json')):
    w = {"id": IID, "q": q, "d": core.degree(D), "D": core.to_str(D),
         "N": [core.to_str(n) for n in N], "Q": core.to_str(Qp)}
    with open(path, 'w') as f:
        json.dump(w, f, indent=1)
    return path


def lin(coeffs):
    """linear form from 6 coefficients"""
    out = {}
    for i, c in enumerate(coeffs):
        if c:
            e = [0]*6; e[i] = 1
            out[tuple(e)] = Fraction(c)
    return out


def var(i):
    e = [0]*6; e[i] = 1
    return {tuple(e): Fraction(1)}


def mono(e):
    return {tuple(e): Fraction(1)}
