import sys, time, json
sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp, served
from qcp import core, P
delta = P(open(qcp.WD + '/delta.txt').read())
Q = qcp.Q
D = core.add(core.add(core.power(Q[0],3),core.power(Q[1],3)),core.power(Q[2],3))
t0=time.time()
res, basis = served.analyse(qcp.F, 1, D, 1, extra={'A': qcp.A, 'delta': delta})
print('join012 odd sector:', res, f'{time.time()-t0:.1f}s', flush=True)
t0=time.time()
res, basis = served.analyse(qcp.F, 1, D, 0)
print('join012 even sector:', res, f'{time.time()-t0:.1f}s', flush=True)
