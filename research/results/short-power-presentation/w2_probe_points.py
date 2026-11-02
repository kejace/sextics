import sys, json, numpy as np
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2_varpro import *
F, A, q, rec = load()
J = json.load(open('kernel_d7_p0.json'))
reps = [core.parse(s) for s in J['basis_rational']]; nums = [[core.parse(c) for c in P] for P in J['numerators_rational']]
rvp = ReducedVarPro(F, A, 7, 0, reps, nums, verbose=True)
pos = {(m['j'], m['m']): k for k, m in enumerate(rvp.fam.members) if m['tag'] == 'J'}
def field_vec(coeffs):   # coeffs: dict (j, m) -> value in unscaled coords
    a = np.zeros(rvp.K)
    for key, v in coeffs.items(): a[pos[key]] += v * rvp.scale[pos[key]]
    return a / np.linalg.norm(a)
def ex(h):   # h * E_X : components (j, h*X_j) for j<3
    return {(j, mon_add(h, tuple(1 if t == j else 0 for t in range(NV)))): 1.0 for j in range(3)}
def ey(h):
    return {(j, mon_add(h, tuple(1 if t == j else 0 for t in range(NV)))): 1.0 for j in range(3, 6)}
tests = {}
for i, name in enumerate(['X0', 'X1', 'X2']):
    h = tuple(1 if t == i else 0 for t in range(NV))
    tests[f'{name}*E_X'] = field_vec(ex(h)); tests[f'{name}*E_Y'] = field_vec(ey(h))
# kernel direction
a = np.zeros(rvp.K); a[0] = 1.0; tests['D7'] = a
for name, a in tests.items():
    rel, f, bb = rvp.rel(a)
    print(f'{name}: f/scale0 = {rel:.3e}  f = {f:.3e}  |b|^2 = {bb:.3e}', flush=True)
# compare with the probe minima if available
import glob
for fjs in glob.glob('rvp_d7_p0_s*_probe.json'):
    C = json.load(open(fjs))['candidates']
    for c in C[:2]:
        a = np.asarray(c['a']); 
        dots = {name: float(abs(a @ v)) for name, v in tests.items()}
        top = sorted(dots.items(), key=lambda kv: -kv[1])[:3]
        # dominant family members
        idx = np.argsort(-np.abs(a))[:5]
        print(f'probe candidate rel={c["rel"]:.2e}: overlaps {top}; dominant members {[(rvp.fam.members[k]["label"], round(float(a[k]),3)) for k in idx]}', flush=True)
