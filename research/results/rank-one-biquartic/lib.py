"""Helper library for the rank-one-biquartic attack: logging + fixed-D runner."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'rank-one-biquartic'
WORK = f'/tmp/hodge/results/{IID}'
LOG = f'{WORK}/log.jsonl'

F, A, q, rec = core.load_instance(IID)
P = core.parse
t = Fraction(5, 7)

f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
L = P('X0**2 + X0*X1 + 2*X1*X2 + 3*X2**2')
M = P('Y0**2 + Y0*Y1 + 2*Y1*Y2 + 3*Y2**2')
Pc = core.mul(core.mul(L, L), M)          # coupling monomial block L^2 M
assert core.sub(F, core.add(core.add(f, g), core.scale(Pc, t))) == {}, 'decomposition check failed'

def eps(Pp, idxs):
    """Radial derivative sum_{i in idxs} Z_i d/dZ_i."""
    out = {}
    for i in idxs:
        out = core.add(out, core.mul({tuple(1 if j == i else 0 for j in range(6)): Fraction(1)}, core.pderiv(Pp, i)))
    return out

RX = eps(F, [0, 1, 2])   # 6f + 4tP
RY = eps(F, [3, 4, 5])   # 6g + 2tP


def log(entry):
    entry.setdefault('ts', time.strftime('%H:%M:%S'))
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(entry, default=str) + '\n')


def run_D(D, expr, family, note='', fast=True, want_exact=True, skip_e1=False, timeout_note=None):
    """Run E1 test and the fixed-D system; log and return result dict."""
    d = core.degree(D)
    t0 = time.time()
    res = {'D': expr, 'd': d, 'family': family, 'note': note}
    if not D:
        res.update(e1=False, consistent=None, rank=0, nullity=0, seconds=0.0, note=note + ' [D=0]')
        log(res); return res
    try:
        if core.parity_of(D) is None:
            res.update(e1=False, consistent=None, rank=0, nullity=0, seconds=0.0, note=note + ' [mixed parity]')
            log(res); return res
        if not skip_e1:
            ok, info = fixed_d.e1_solvable(F, A, q, D)
            res['e1'] = bool(ok)
            res['e1_info'] = info
            if not ok:
                res.update(consistent=False, rank=info['rank'], nullity=info['cols'] - info['rank'],
                           seconds=round(time.time() - t0, 2))
                res['note'] = note + ' [E1 fails]'
                log(res); return res
        else:
            res['e1'] = None
        if fast:
            out = fixed_d.solve_fixed_D_fast(F, A, q, D, want_exact=want_exact)
            res.update(consistent=bool(out['consistent']), rank=int(out['rank']), nullity=int(out['nullity']),
                       rows=out['rows'], cols=out['cols'])
            if out.get('consistent'):
                res['exact'] = out.get('exact')
                if out.get('exact'):
                    res['N'] = [core.to_str(n) for n in out['N']]
                    res['Q'] = core.to_str(out['Q'])
        else:
            status, data = fixed_d.solve_fixed_D(F, A, q, D)
            res['status'] = status
            if status == 'inconsistent':
                res.update(consistent=False, rank=data['rank'], nullity=data['nullity'], rows=data['rows'], cols=data['cols'])
            elif status == 'solved':
                res.update(consistent=True, rank=data['rank'], nullity=data['nullity'], exact=True,
                           N=[core.to_str(n) for n in data['N']], Q=core.to_str(data['Q']))
            else:
                res.update(consistent=True, rank=-1, nullity=-1, exact=False)
    except Exception as ex:
        res.update(consistent=None, rank=-1, nullity=-1, note=note + f' [error: {ex!r}]')
    res['seconds'] = round(time.time() - t0, 2)
    log(res)
    if res.get('consistent') and res.get('exact'):
        with open(f'{WORK}/witness-candidate.json', 'w') as fh:
            json.dump({'instance': IID, 'd': d, 'D': core.to_str(D), 'N': res['N'], 'Q': res['Q']}, fh, indent=1)
        print('!!!!!!!!!! WITNESS FOUND !!!!!!!!!!', expr, file=sys.stderr)
    return res


def show(res):
    keys = ['D', 'd', 'family', 'e1', 'consistent', 'rank', 'nullity', 'rows', 'cols', 'seconds', 'note']
    print(' | '.join(f'{k}={res.get(k)}' for k in keys if k in res), flush=True)
