"""E5: degree-6 Hessian relaxation on the FULL E1-admissible space (all sectors, all parities).
Quadrics q_{psi,s}(c) = coordinates of Phi(sum c_k D_k) in R_17 (standard monomials per sector).
Writes a Singular input (char p) and an msolve input; the variety must be computed by a GB run."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
d = int(sys.argv[1]) if len(sys.argv) > 1 else 6
p = P2
t0 = time.time()
JS = JacobianSectors(F2, p)
Ds = []; keys = []
for key in all_sector_keys(JS.gens, d):
    r = JS.e1_kernel(A, d, key)
    for Dk in r['basis']:
        Ds.append(Dk); keys.append(key)
print('E1 space dim', len(Ds), f'{time.time()-t0:.0f}s', flush=True)
N0, B = relaxation_data(JS, A, Ds)
print('bilinear forms done', f'{time.time()-t0:.0f}s', flush=True)
# reduce each B(k,l) to normal-form coordinates in R_{2d+... } degree = d + (d+7) - 1 = 2d+6
degB = 2 * d + 6
kA = key_of_poly(A, JS.gens)
quadrics = {}   # (psi, s) -> {(k,l): coeff}
nz = 0
for (k, l), Bkl in B.items():
    if not Bkl:
        continue
    psi = key_of_poly(Bkl, JS.gens)
    C = JS.normal_form_coords([Bkl], degB, psi)
    for s in range(C.ncols()):
        v = int(C[0, s])
        if v:
            quadrics.setdefault((psi, s), {})[(k, l)] = v; nz += 1
print('nonzero quadric coefficients', nz, 'quadrics', len(quadrics), f'{time.time()-t0:.0f}s', flush=True)
# sanity: F2 is in the E1 space and Phi(F2) must vanish: express F2 in the basis -> need coordinates; instead verify directly
# Phi(F2) via relaxation_data on [F2]
_, BF = relaxation_data(JS, A, [F2])
CF = JS.normal_form_coords([BF[(0, 0)]], degB, key_of_poly(BF[(0, 0)], JS.gens)) if BF[(0, 0)] else None
phiF_zero = (CF is None) or all(int(CF[0, s]) == 0 for s in range(CF.ncols()))
print('Phi(F) == 0 in R:', phiF_zero, flush=True)
# write ideal files
vars_ = [f'c{k}' for k in range(len(Ds))]
polys = []
for (psi, s), terms in quadrics.items():
    mon = []
    for (k, l), v in terms.items():
        coef = v if k == l else (2 * v) % p
        mon.append(f'{coef}*c{k}*c{l}' if k != l else f'{coef}*c{k}^2')
    polys.append('+'.join(mon))
base = f'/tmp/hodge/results/cycle-q2/wave2-excl/receipts/hess{d}_p{p}'
with open(base + '.sing', 'w') as f:
    f.write(f'ring r = {p},({",".join(vars_)}),dp;\n')
    f.write('ideal I = ' + ',\n'.join(polys) + ';\n')
    f.write('option(redSB);\nideal G = std(I);\nprint("size GB=" + string(size(G)));\nprint("dim=" + string(dim(G)));\nprint("degree=" + string(mult(G)));\n')
    f.write('print("hilbPoly:"); hilbPoly(G);\nprint("GB:"); print(G);\nquit;\n')
with open(base + '.ms', 'w') as f:
    f.write(','.join(vars_) + '\n' + str(p) + '\n')
    f.write(',\n'.join(pl.replace('^', '^') for pl in polys) + '\n')
json.dump(dict(d=d, p=p, basis=[to_str(Dk) for Dk in Ds], keys=[list(k) for k in keys], n_quadrics=len(quadrics), phiF_zero=phiF_zero),
          open(base + '_basis.json', 'w'))
rec = dict(exp=f'E5-hessian-relaxation-setup', engine=ENGINE_TAG, d=d, p=p, e1_dim=len(Ds), n_quadrics=len(quadrics), nonzero_coeffs=nz,
           phiF_zero=phiF_zero, files=[base + '.sing', base + '.ms'], seconds=round(time.time() - t0, 1))
print(json.dumps(rec)); xlib.log(rec)
