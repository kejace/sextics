"""Helpers for the symmetric-sextic-residual attack: logging + families."""
import sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'symmetric-sextic-residual'
RES = f'/tmp/hodge/results/{IID}'
LOG = f'{RES}/log.jsonl'
P = core.parse

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def test_D(F, A, q, D, family, note='', full=True, fast=False, verbose=False, label=None):
    """E1 test, then (optionally) full fixed-D system. Logs and returns dict."""
    expr = label or core.to_str(D)
    d = core.degree(D)
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    out = {'D': expr, 'd': d, 'e1': ok, 'e1_info': info}
    if not ok or not full:
        log(expr, d, family, ok, None, info['rank'], info['cols'] - info['rank'], t1 - t0,
            note + f" | E1 only: rank={info['rank']} rank_aug={info['rank_aug']} rows={info['rows']} cols={info['cols']}")
        out['consistent'] = None
        return out
    if fast:
        r = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
        t2 = time.time()
        out.update(consistent=r['consistent'], rank=r['rank'], nullity=r['nullity'], rows=r['rows'], cols=r['cols'])
        if r['consistent'] and r.get('exact'):
            out['N'], out['Q'] = r['N'], r['Q']
        log(expr, d, family, ok, r['consistent'], r['rank'], r['nullity'], t2 - t0,
            note + f" | fast full system rows={r['rows']} cols={r['cols']} exact={r.get('exact')}")
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
        t2 = time.time()
        cons = status != 'inconsistent'
        out.update(status=status, consistent=cons, rank=data.get('rank'), nullity=data.get('nullity'),
                   rows=data.get('rows'), cols=data.get('cols'))
        if status == 'solved':
            out['N'], out['Q'] = data['N'], data['Q']
        log(expr, d, family, ok, cons, data.get('rank'), data.get('nullity'), t2 - t0,
            note + f" | full system status={status} rows={data.get('rows')} cols={data.get('cols')}")
    return out

# --- symmetric-function helpers ---------------------------------------------
def var(i):
    e = [0]*6; e[i] = 1
    return {tuple(e): Fraction(1)}

X = [var(0), var(1), var(2)]
Y = [var(3), var(4), var(5)]

def e_k(vs, k):
    out = {}
    for c in itertools.combinations(vs, k):
        t = {(0,)*6: Fraction(1)}
        for v in c:
            t = core.mul(t, v)
        out = core.add(out, t)
    return out

def p_k(vs, k):
    out = {}
    for v in vs:
        out = core.add(out, core.power(v, k))
    return out

def lin(coeffs):
    """linear form sum c_i Z_i"""
    out = {}
    for i, c in enumerate(coeffs):
        if c:
            out = core.add(out, core.scale(var(i), Fraction(c)))
    return out

def permute(Pol, perm):
    """apply variable permutation: variable i -> perm[i]"""
    out = {}
    for e, c in Pol.items():
        e2 = [0]*6
        for i in range(6):
            e2[perm[i]] += e[i]
        out[tuple(e2)] = out.get(tuple(e2), 0) + c
    return {e: c for e, c in out.items() if c}

def signs(Pol, sg):
    out = {}
    for e, c in Pol.items():
        s = 1
        for i in range(6):
            if sg[i] < 0 and e[i] % 2: s = -s
        out[e] = c * s
    return out
