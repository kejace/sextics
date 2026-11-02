import sys, time; sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import w2lib as W, numpy as np
from hodge import core
F, A, q = W.load()
G = W.group('sgnX')
fam = W.family(F, A, q, 7, 0, 'sgnX')
basis = fam['kernel'] + fam['jac']
rng = np.random.default_rng(3); a = rng.standard_normal(len(basis))
t0 = time.time(); c1 = W.CharVarProSlow(F, A, q, basis, G); t1 = time.time() - t0
t0 = time.time(); c2 = W.CharVarPro(F, A, q, basis, G); t2 = time.time() - t0
print('cvp slow/fast build', round(t1,1), round(t2,1), 'values', c1.value(a)[0], c2.value(a)[0], 'Mk diff', max(abs((c1.Mk[k]-c2.Mk[k])).max() for k in range(len(basis))), 'V diff', np.abs(c1.V-c2.V).max(), 'MQ diff', np.abs(c1.MQ-c2.MQ).max())
N0 = W.e1_lift_char(F, A, q, fam['kernel'][0], G)
pairs = [(fam['kernel'][0], N0)] + W.jac_family_with_lifts(F, A, q, 7, 0, G)
t0 = time.time(); k1 = W.KoszulVarProSlow(F, A, q, pairs, G); t1 = time.time() - t0
t0 = time.time(); k2 = W.KoszulVarPro(F, A, q, pairs, G); t2 = time.time() - t0
print('kvp slow/fast build', round(t1,1), round(t2,1), 'values', k1.value(a)[0], k2.value(a)[0], 'V diff', np.abs(k1.V-k2.V).max(), 'nC', k1.nC, k2.nC, 'MQ diff', np.abs(k1.MQ-k2.MQ).max())
# float lifts equivalence
N0f, rr = W.e1_lift_float(F, A, q, fam['kernel'][0], G)
pairsf = [(fam['kernel'][0], N0f)] + W.jac_family_with_lifts(F, A, q, 7, 0, G)
k3 = W.KoszulVarPro(F, A, q, pairsf, G)
print('float lift relres', rr, 'value with float lift', k3.value(a)[0])
# family_float at d=9 sgnX vs exact (K=19)
polys, info = W.family_float(F, A, q, 9, 0, 'sgnX', K_expected=19, verbose=True)
fam9 = W.family(F, A, q, 9, 0, 'sgnX')
# check: each exact basis element lies in the float span (project)
mons = {}
for P in polys + fam9['kernel'] + fam9['jac']:
    for e in P: mons.setdefault(e, len(mons))
Mf = np.zeros((len(mons), len(polys)))
for j, P in enumerate(polys):
    for e, v in P.items(): Mf[mons[e], j] = v
Qm, _ = np.linalg.qr(Mf)
worst = 0
for P in fam9['kernel'] + fam9['jac']:
    v = np.zeros(len(mons))
    for e, c in P.items(): v[mons[e]] = float(c)
    res = np.linalg.norm(v - Qm @ (Qm.T @ v)) / np.linalg.norm(v); worst = max(worst, res)
print('exact family inside float span: worst rel residual', worst)
# timing at d=9 sgnX kvp
pairs9 = []
for P in fam9['kernel']:
    N0, rr = W.e1_lift_float(F, A, q, P, G); pairs9.append((P, N0))
pairs9 += W.jac_family_with_lifts(F, A, q, 9, 0, G)
t0 = time.time(); kv = W.KoszulVarPro(F, A, q, pairs9, G, verbose=True); 
a9 = rng.standard_normal(kv.K)
t0 = time.time(); kv.value(a9); print('kvp d9 eval', round(time.time()-t0, 3))
t0 = time.time(); kv.value_and_grad(kv.P.T @ a9); print('kvp d9 eval+grad', round(time.time()-t0, 3))
t0 = time.time(); cv = W.CharVarPro(F, A, q, [P for P, _ in pairs9], G, verbose=True)
t0 = time.time(); cv.value(a9); print('cvp d9 eval', round(time.time()-t0, 3))
