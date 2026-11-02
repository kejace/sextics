import sys, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from hlib import *
F, A, q, rec = load()
print('F terms', len(F), 'A', core.to_str(A), 'q', q, 'parity A', core.parity_of(A), 'parity F', core.parity_of(F))
# sanity: F == eX + eY + c((L+M)^6 + (L-M)^6) ?
Fp = core.add(core.add(eX, eY), core.scale(core.add(core.power(core.add(L, M), 6), core.power(core.sub(L, M), 6)), c))
print('F == eX+eY+c((L+M)^6+(L-M)^6):', core.sub(F, Fp) == {})
# (0) D = 1
status, data = run(F, A, q, P('1'), 'D=1', note='is A in J(F)?')
print('D=1:', status, data)
