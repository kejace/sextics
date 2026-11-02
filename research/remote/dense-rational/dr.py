"""Helpers for the dense-rational attack. Not part of the repo."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import NV, monomials, shift, mul, add, sub, scale, pderiv, power, yparity
from fractions import Fraction
import flint
from flint import nmod_mat

IID = 'dense-rational'
LOG = '/tmp/hodge/results/dense-rational/log.jsonl'
P = core.PRIMES[0]

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def test_D(F, A, q, D, family, name=None, full=True, fast=True, note=''):
    """E1 test, then (optionally) full fixed-D system.  Logs result.  Returns dict."""
    name = name or core.to_str(D)
    d = core.degree(D)
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    res = {'name': name, 'd': d, 'e1': ok, 'e1_info': info, 't_e1': round(t1 - t0, 1)}
    if not ok or not full:
        log(name, d, family, ok, None if ok else False, info['rank'], info['cols'] - info['rank'], t1 - t0,
            note=(note + ' E1 only' if ok else note + ' E1 fails (rank %d, rank_aug %d)' % (info['rank'], info['rank_aug'])).strip())
        return res
    t2 = time.time()
    if fast:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True)
        cons, rank, nullity = out['consistent'], out['rank'], out['nullity']
        res.update(out)
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True)
        cons = status == 'solved'
        rank, nullity = data.get('rank'), data.get('nullity')
        res['status'] = status; res['data'] = data
    t3 = time.time()
    res['t_full'] = round(t3 - t2, 1)
    log(name, d, family, True, cons, rank, nullity, t3 - t0, note=note + (' CONSISTENT!!' if cons else ''))
    return res

def psym(P_):
    """Signed representation of a mod-p poly for printing."""
    return {k: (v if v < P // 2 else v - P) for k, v in P_.items()}

def hilbert_counts():
    """dim S_d^{+/-}, dim R_d^{+/-} for the Jacobian ring of a smooth sextic (parity-graded)."""
    import numpy as np
    # H(t) = (1+t+..+t^4)^6 ; Theta(t) = (1+t^2+..+t^8)^3
    def polymul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    H = [1]
    for _ in range(6): H = polymul(H, [1, 1, 1, 1, 1])
    Th = [1]
    for _ in range(3): Th = polymul(Th, [1, 0, 1, 0, 1, 0, 1, 0, 1])
    H += [0] * 40; Th += [0] * 40
    R = {}
    for d in range(0, 30):
        R[d] = ((H[d] + Th[d]) // 2, (H[d] - Th[d]) // 2)
    S = {}
    for d in range(0, 40):
        S[d] = (len(monomials(d, 0)), len(monomials(d, 1)))
    return S, R
