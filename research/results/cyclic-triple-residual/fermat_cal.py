"""Fermat calibration (t=0) with the sector solver: which D serve X0^3 Y0^3 on the Fermat sextic?"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
def log(e):
    with open(f'{WD}/log.jsonl', 'a') as f: f.write(json.dumps(e) + '\n')
F0 = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
A = core.parse('X0**3*Y0**3')
cases = [
 ('X0**6 + X1**6 + X2**6', 'fermat:join3-X', 'join of plane sextic curves C_X * C_Y: vertex P^2_Y, base genus-10 curve'),
 ('Y0**6 + Y1**6 + Y2**6', 'fermat:join3-Y', 'same divisor (F - D)'),
 ('X0', 'fermat:reflection-hyperplane', 'branch divisor of X0 -> -X0'),
 ('X0**2', 'fermat:reflection-hyperplane-sq', ''),
 ('X0**3', 'fermat:reflection-hyperplane-cube', ''),
 ('X0*Y0', 'fermat:two-hyperplanes', ''),
 ('X0**3*Y0**3', 'fermat:A-itself', ''),
 ('X0**6 + Y0**6', 'fermat:cone-hyperplanes', 'crosscheck of the engine-2 result (inconsistent)'),
 ('X0**6 + Y0**6 + Y1**6', 'fermat:sum3-mixed', ''),
 ('X1**6 + X2**6', 'fermat:sum2', ''),
 ('X1**6 + X2**6 + Y1**6 + Y2**6', 'fermat:sum4-complement', ''),
 ('X0**6 + X1**6 + X2**6 + Y0**6', 'fermat:cone4-sextic', 'crosscheck (inconsistent)'),
 ('X0**3 + X1**3 + X2**3', 'fermat:cone3-cubic', 'crosscheck (E1 fails)'),
 ('X0**2 + X1**2 + X2**2', 'fermat:cone3-quadric', ''),
 ('X0**6 + X1**6 + X2**6 + 2*Y0**6 + 3*Y1**6 + 5*Y2**6', 'fermat:generic-sum6', ''),
]
for Dexpr, fam, note in cases:
    D = core.parse(Dexpr)
    t0 = time.time()
    out = test(F0, A, 1, D, exact=True, verbose=False)
    ent = {'engine': 2, 'D': Dexpr, 'd': out['d'], 'family': fam, 'e1': None, 'consistent': out['consistent'], 'rank': out['rank'],
           'nullity': out['nullity'], 'seconds': round(time.time() - t0, 2), 'tag': 'fermat',
           'note': 'FERMAT-CALIBRATION (sector solver, |G|=%d, rows=%d cols=%d rank_aug=%d, exact over Q) %s' % (out['G_order'], out['rows'], out['cols'], out['rank_aug'], note)
                   + (' | SOLVED on Fermat, verified=%s' % out.get('verified') if out['consistent'] else '')}
    if out.get('verified'): ent['witness'] = {'D': Dexpr, 'N': out['N'], 'Q': out['Q']}
    log(ent)
    print(f"d={out['d']:2d} cons={str(out['consistent']):5s} rank={out['rank']:5d}/{out['rank_aug']:5d} cols={out['cols']:5d} rows={out['rows']:5d} {ent['seconds']:6.1f}s verified={out.get('verified')} {fam:34s} {Dexpr}", flush=True)
    if out.get('verified'):
        print('   N =', out['N']); print('   Q =', out['Q'])
