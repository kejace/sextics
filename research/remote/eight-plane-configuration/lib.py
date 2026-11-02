"""Helper library for the eight-plane-configuration attack."""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'eight-plane-configuration'
WD = f'/tmp/hodge/results/{IID}'
LOG = f'{WD}/log.jsonl'

F, A, q, rec = core.load_instance(IID)

def log(entry):
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')

def test_D(D, name, family, note='', fast=True, want_exact=True, e1_only=False, verbose=False, p=None):
    """Run E1 test then the full fixed-D system. Logs and returns a dict."""
    if isinstance(D, str):
        Dexpr = D; D = core.parse(D)
    else:
        Dexpr = name
    d = core.degree(D)
    t0 = time.time()
    entry = {'D': Dexpr if len(Dexpr) < 400 else name, 'name': name, 'd': d, 'family': family, 'e1': None,
             'consistent': None, 'rank': None, 'nullity': None, 'seconds': None, 'note': note}
    try:
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
        entry['e1'] = bool(ok)
        entry['e1_info'] = {k: int(v) for k, v in info.items()}
        if ok and not e1_only:
            if fast:
                out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=verbose, want_exact=want_exact, p=p)
                entry['consistent'] = bool(out['consistent']); entry['rank'] = int(out['rank']); entry['nullity'] = int(out['nullity'])
                entry['rows'] = int(out['rows']); entry['cols'] = int(out['cols'])
                if out['consistent']:
                    entry['exact'] = out.get('exact')
                    if out.get('exact'):
                        entry['note'] += ' SOLVED-EXACT'
                        save_witness(D, out['N'], out['Q'], name)
            else:
                status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=verbose)
                entry['status'] = status
                entry['rank'] = int(data.get('rank', -1)); entry['nullity'] = int(data.get('nullity', -1))
                entry['consistent'] = status != 'inconsistent'
                if status == 'solved':
                    entry['note'] += ' SOLVED-EXACT'
                    save_witness(D, data['N'], data['Q'], name)
    except Exception as ex:
        entry['note'] += f' ERROR {type(ex).__name__}: {ex}'
    entry['seconds'] = round(time.time() - t0, 2)
    log(entry)
    return entry

def save_witness(D, N, Q, name):
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    w = {'instance': IID, 'name': name, 'verified': bool(ok), 'info': str(info), 'd': core.degree(D),
         'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}
    with open(f'{WD}/witness-candidate.json', 'w') as f:
        json.dump(w, f, indent=1)
    print('!!!! WITNESS SAVED', name, ok, info, file=sys.stderr)

def P(s):
    return core.parse(s)

X0 = P('X0'); X1 = P('X1'); X2 = P('X2'); Y0 = P('Y0'); Y1 = P('Y1'); Y2 = P('Y2')
Q0 = P('X0**2+Y0**2'); Q1 = P('X1**2+Y1**2'); Q2 = P('X2**2+Y2**2')
