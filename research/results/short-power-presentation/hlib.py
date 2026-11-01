"""Helpers for the short-power-presentation attack: logging + candidate denominators."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'short-power-presentation'
RES = '/tmp/hodge/results/short-power-presentation'
LOG = os.path.join(RES, 'log.jsonl')

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
           "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec

def load():
    F, A, q, rec = core.load_instance(IID)
    return F, A, q, rec

# linear forms
X0 = core.parse('X0'); X1 = core.parse('X1'); X2 = core.parse('X2')
Y0 = core.parse('Y0'); Y1 = core.parse('Y1'); Y2 = core.parse('Y2')
L = core.parse('X0 + 2*X1 + 3*X2')
M = core.parse('Y0 + 3*Y1 + 5*Y2')
eX = core.parse('X0**6 + X1**6 + X2**6')
eY = core.parse('Y0**6 + Y1**6 + Y2**6')
c = Fraction(5, 7**7)

def P(s):
    return core.parse(s)

def run(F, A, q, D, family, note='', e1_only=False, fast=False, verbose=True):
    """E1 test, then full fixed-D solve. Returns (status, data)."""
    d = core.degree(D)
    Dstr = core.to_str(D)
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    if verbose:
        print(f'[E1] d={d} {family}: e1={ok} {info} ({t1-t0:.1f}s)', file=sys.stderr, flush=True)
    if not ok:
        log(Dstr, d, family, False, False, info['rank'], info['cols'] - info['rank'], t1 - t0, note + ' E1 fails')
        return 'e1-fail', info
    if e1_only:
        log(Dstr, d, family, True, None, info['rank'], info['cols'] - info['rank'], t1 - t0, note + ' E1 only')
        return 'e1-ok', info
    t2 = time.time()
    if fast:
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose)
        t3 = time.time()
        status = 'solved' if out.get('exact') else ('inconsistent' if not out['consistent'] else 'lift-failed')
        log(Dstr, d, family, True, out['consistent'], out['rank'], out['nullity'], t3 - t0, note + f" fast rows={out['rows']} cols={out['cols']}")
        return status, out
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
    t3 = time.time()
    cons = status != 'inconsistent'
    log(Dstr, d, family, True, cons, data.get('rank'), data.get('nullity'), t3 - t0,
        note + f" status={status} rows={data.get('rows')} cols={data.get('cols')}")
    if verbose:
        print(f'[solve] {family}: {status} ({t3-t2:.1f}s)', file=sys.stderr, flush=True)
    return status, data

def served(F, q, D, verbose=False):
    """target_image wrapper: (dim, basis, S)"""
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D, verbose=verbose)
    return dim, basis, S, time.time() - t0
