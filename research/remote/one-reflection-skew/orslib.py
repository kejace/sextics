"""Helper library for one-reflection-skew experiments: logging + wrappers."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

INST = 'one-reflection-skew'
WD = f'/tmp/hodge/results/{INST}'
LOG = f'{WD}/log.jsonl'
F, A, q, rec = core.load_instance(INST)

def log(D_expr, d, family, e1, consistent, rank, nullity, seconds, note=''):
    rec_ = {"D": D_expr, "d": d, "family": family, "e1": e1, "consistent": consistent,
            "rank": rank, "nullity": nullity, "seconds": round(seconds, 2), "note": note, "engine": 2}
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec_) + '\n')
    return rec_

def test(D_expr, family, note='', full=True, fast=False, want_exact=True, Dpoly=None, verbose=False):
    """E1 test, then (optionally) full fixed-D solve.  Appends to log.  Returns dict."""
    D = Dpoly if Dpoly is not None else core.parse(D_expr)
    d = core.degree(D)
    t0 = time.time()
    e1, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    out = {'D': D_expr, 'd': d, 'e1': e1, 'e1_info': info, 'consistent': None, 'rank': None,
           'nullity': None, 'status': None}
    if not e1 or not full:
        out['seconds'] = t1 - t0
        log(D_expr, d, family, e1, None, info['rank'], info['cols'] - info['rank'], t1 - t0,
            note + f' | E1 rank={info["rank"]} aug={info["rank_aug"]} rows={info["rows"]} cols={info["cols"]}')
        return out
    if fast:
        r = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=want_exact)
        out.update({'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
                    'rows': r['rows'], 'cols': r['cols']})
        if r['consistent'] and r.get('exact'):
            out['status'] = 'solved'; out['N'] = r['N']; out['Q'] = r['Q']
        elif r['consistent']:
            out['status'] = 'consistent-modp'
        else:
            out['status'] = 'inconsistent'
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
        out['status'] = status
        out.update({'consistent': status != 'inconsistent', 'rank': data.get('rank'),
                    'nullity': data.get('nullity'), 'rows': data.get('rows'), 'cols': data.get('cols')})
        if status == 'solved':
            out['N'] = data['N']; out['Q'] = data['Q']
    t2 = time.time()
    out['seconds'] = t2 - t0
    log(D_expr, d, family, e1, out['consistent'], out['rank'] or 0, out['nullity'] or 0, t2 - t0,
        note + f' | status={out["status"]}')
    if out['status'] == 'solved':
        print('!!!!!!!! WITNESS FOUND !!!!!!!!', D_expr, file=sys.stderr)
        save_witness(D, out['N'], out['Q'])
    return out

def save_witness(D, N, Q):
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    w = {'instance': INST, 'q': q, 'd': core.degree(D), 'D': core.to_str(D),
         'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q), 'verified': ok, 'info': str(info)}
    with open(f'{WD}/witness-candidate.json', 'w') as f:
        json.dump(w, f, indent=1)
    return ok

# convenience polynomial constructors
X = [core.parse(v) for v in core.VARS]
def P(s): return core.parse(s)
def Fi(i): return core.pderiv(F, i)
