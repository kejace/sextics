import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import c, L, M, eX, eY, ppow
F, A, q, rec = core.load_instance(IID)
L2, M2 = core.mul(L, L), core.mul(M, M)
# D_res = e(X)(L^4 - L^2 M^2 + M^4) + 2c L^6 (L^4 + 14 L^2 M^2 + M^4)
P1 = core.sub(core.add(core.mul(L2, L2), core.mul(M2, M2)), core.mul(L2, M2))
P2 = core.add(core.add(core.mul(L2, L2), core.scale(core.mul(L2, M2), 14)), core.mul(M2, M2))
Dres = core.add(core.mul(eX, P1), core.scale(core.mul(ppow(L, 6), P2), 2 * c))
Dres_int = core.scale(Dres, 7**7)   # integral
print('deg', core.degree(Dres), 'terms', len(Dres), 'parity', core.parity_of(Dres))
# binary-form reduction: 10c L^2 M^2 (L^2+M^2)^2 (L^2-M^2)
Bf = core.scale(core.mul(core.mul(core.mul(L2, M2), core.power(core.add(L2, M2), 2)), core.sub(L2, M2)), 10 * c)
for p in core.PRIMES[:3]:
    vecs, info = e1_kernel(F, A, [('Dres', Dres_int)], 10, 0, p)
    jv, infoj = in_J_family(F, [('Dres', Dres_int)], 10, 0, p)
    # is Dres - Bf in J_10 ? (checks the collapse e(X) = -L h_L/6 mod J on this example)
    jv2, _ = in_J_family(F, [('diff', core.scale(core.sub(Dres, Bf), 7**7))], 10, 0, p)
    vb, infob = e1_kernel(F, A, [('Bf', core.scale(Bf, 7**7))], 10, 0, p)
    r = {'what': 'E1-Dres', 'd': 10, 'parD': 0, 'prime': p, 'E1_Dres': info['kernel_dim'], 'Dres_in_J10': infoj['kernel_dim'],
         'Dres_minus_binaryform_in_J10': jv2 and len(jv2) or 0, 'E1_binaryform': infob['kernel_dim'], 'hilbert_ok': info['hilbert_ok'] and infoj['hilbert_ok'],
         'result': 'A*Dres in J: %s ; Dres in J_10: %s ; Dres == 10c L^2M^2(L^2+M^2)^2(L^2-M^2) mod J_10: %s' % (bool(info['kernel_dim']), bool(infoj['kernel_dim']), bool(jv2)),
         'file': 'w4_dres.py'}
    log(r); print(json.dumps(r))
