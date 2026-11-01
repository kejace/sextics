"""Harness for cyclic-triple-residual experiments (exact, mod-p via fixed_d)."""
import sys, time, json, os, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction

IID = 'cyclic-triple-residual'
WD = '/tmp/hodge/results/cyclic-triple-residual'
LOG = f'{WD}/log.jsonl'
F, A, q, rec = core.load_instance(IID)
V = core.VARS
t = Fraction(5, 7)

def P(s):
    return core.parse(s)

def log(entry):
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry) + '\n')

def test_D(D, family, note='', Dexpr=None, F_=None, A_=None, q_=None, mode='auto', second_prime=True, verbose=False, tag=None):
    """Run E1 test then full fixed-D solve. Returns dict. Logs to log.jsonl."""
    F_ = F if F_ is None else F_; A_ = A if A_ is None else A_; q_ = q if q_ is None else q_
    if isinstance(D, str):
        Dexpr = D; D = P(D)
    if Dexpr is None:
        Dexpr = core.to_str(D)
    d = core.degree(D)
    if not core.is_homogeneous(D) or core.parity_of(D) is None:
        ent = {'engine': 2, 'D': Dexpr[:400], 'd': d, 'family': family, 'e1': None, 'consistent': None, 'rank': -1, 'nullity': -1,
               'seconds': 0.0, 'note': note + ' SKIPPED: not homogeneous or mixed parity'}
        log(ent); return ent
    t0 = time.time()
    e1, info1 = fixed_d.e1_solvable(F_, A_, q_, D)
    res = {'engine': 2, 'D': Dexpr[:400], 'd': d, 'family': family, 'e1': bool(e1), 'consistent': None, 'rank': int(info1['rank']),
           'nullity': int(info1['cols'] - info1['rank']), 'seconds': 0.0, 'note': note, 'e1_rows': info1['rows'], 'e1_cols': info1['cols']}
    if tag: res['tag'] = tag
    if not e1:
        res['seconds'] = round(time.time() - t0, 2)
        res['note'] = (note + ' | E1 fails (A*D not in J(F)) rank=%d rank_aug=%d' % (info1['rank'], info1['rank_aug'])).strip(' |')
        log(res); return res
    # E1 passes: full system
    if mode == 'auto':
        mode = 'full' if d <= 4 else 'fast'
    if mode == 'full':
        status, data = fixed_d.solve_fixed_D(F_, A_, q_, D, verbose=verbose, max_primes=6)
        res['consistent'] = status != 'inconsistent'
        res['rank'] = int(data.get('rank', -1)); res['nullity'] = int(data.get('nullity', -1))
        res['status'] = status
        if status == 'inconsistent' and second_prime and d <= 4:
            status2, data2 = fixed_d.solve_fixed_D(F_, A_, q_, D, primes=core.PRIMES[1:], verbose=verbose, max_primes=1)
            res['status2'] = status2
            res['consistent'] = status2 != 'inconsistent'
        if status == 'solved':
            res['solved'] = True
            N, Q = data['N'], data['Q']
            ok, vinfo = core.verify_witness(F_, A_, q_, D, N, Q)
            res['verified'] = bool(ok)
            res['witness'] = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}
    else:
        out = fixed_d.solve_fixed_D_fast(F_, A_, q_, D, verbose=verbose, want_exact=True)
        res['consistent'] = bool(out['consistent']); res['rank'] = int(out['rank']); res['nullity'] = int(out['nullity'])
        res['rows'] = out['rows']; res['cols'] = out['cols']
        if (not out['consistent']) and second_prime:
            out2 = fixed_d.solve_fixed_D_fast(F_, A_, q_, D, p=core.PRIMES[1], verbose=verbose, want_exact=False)
            res['consistent2'] = bool(out2['consistent']); res['rank2'] = int(out2['rank'])
        if out['consistent'] and out.get('exact'):
            res['solved'] = True
            N, Q = out['N'], out['Q']
            ok, vinfo = core.verify_witness(F_, A_, q_, D, N, Q)
            res['verified'] = bool(ok)
            res['witness'] = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)}
        elif out['consistent']:
            res['note'] += ' | CONSISTENT mod p but exact lift failed'
    res['seconds'] = round(time.time() - t0, 2)
    if res.get('solved'):
        res['note'] += ' | *** SOLVED ***'
        with open(f'{WD}/witness-candidate-{int(time.time())}.json', 'w') as f:
            json.dump({'instanceId': IID, 'family': family, 'verified_by_core': res.get('verified'), **res['witness']}, f, indent=1)
        print('\n*** SOLVED ***', Dexpr, 'verified=', res.get('verified'), '\n', flush=True)
    log(res)
    return res

def fmt(r):
    return f"d={r['d']:2d} e1={str(r['e1']):5s} cons={str(r.get('consistent')):5s} rank={r['rank']:5d} null={r['nullity']:4d} {r['seconds']:7.1f}s  {r['family']:22s} {r['D'][:70]}"

# ---------- structural pieces ----------
def pieces():
    X0,X1,X2,Y0,Y1,Y2 = [P(v) for v in V]
    C = P('X0**3 + X0*Y1*Y2 + X1*Y0*Y1')
    H = P('X0**6 + X1**6 + Y0**6 + Y1**6 + Y2**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0**2*Y0*Y1*Y2**2/7 + 5*X0*X1*Y1**3*Y2/7')
    Z = P('X2**3')
    U = core.add(Z, core.scale(C, t/2))
    Hp = core.sub(H, core.scale(core.mul(C, C), t*t/4))
    Delta = core.sub(core.scale(core.mul(C, C), t*t), core.scale(H, 4))
    return dict(C=C, H=H, Z=Z, U=U, Hp=Hp, Delta=Delta)

if __name__ == '__main__':
    pass
