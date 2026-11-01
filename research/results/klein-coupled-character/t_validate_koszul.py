import sys, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc2 import *
for iid in ['cone-ladder-linear', 'two-reflection-cubic-surface-d2', 'f3-block-x0-quartic-strict', 'hesse-radial-hesse-cofactor-d4-rational']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    t = time.time()
    e = full_test(F, A, D, 'witness '+iid, 'validate', weights=False, want_log=False)
    print(iid, 'd=', d, 'e1', e['e1'], 'consistent', e['consistent'], 'rank', e['rank'], 'nullity', e['nullity'], f"{time.time()-t:.1f}s")
    if e['consistent']:
        r = lift_certificate(F, A, D, weights=False)
        print('   lifted & verified:', None if r is None else (r['ok'], r['info']))
# our F: D = F_0 and Jhat
F, A, q = load()
Fi = [core.pderiv(F, i) for i in range(6)]
for i in range(3):
    e = full_test(F, A, Fi[i], f'F_{i}', 'polar', note='partial derivative dF/dX%d' % i, verbose=True)
    print('F_%d' % i, e['consistent'], e['rank'], e['nullity'], e['seconds'])
Jhat = core.parse('7*X0*X2**2*(9*(Y0**6 - Y1**6 + Y2**6) - 9*(' + open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0] + ')) + 15*X1*Y2*Y0*(2*Y0**6 - 4*Y1**6 + 8*Y2**6 - 2*(' + open('/home/kejace/devel/sextics/data/klein-coupled-character/instance.rec').read().split('F: ')[1].split('\n')[0] + '))')
print('Jhat deg', core.degree(Jhat), 'weight', sigma_weight(Jhat), 'terms', len(Jhat))
e = full_test(F, A, Jhat, 'Jhat (fqlx)', 'companion', note='fqlx degree-9 denominator, character 2', verbose=True)
print('Jhat', e['consistent'], e['rank'], e['nullity'], e['rows'], e['cols'], e['seconds'])
