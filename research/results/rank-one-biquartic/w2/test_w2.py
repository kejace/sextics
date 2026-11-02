import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/w2'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np
import w2lib as W
from w2lib import *
from hodge import core
P = core.parse
t = Fraction(5, 7)
f = P('X0**6 + 5*X0**5*X1/7 + 5*X0**2*X1**3*X2/7 + 15*X0*X2**5/7 + X1**6 + 10*X1**5*X2/7 + X2**6')
g = P('Y0**6 + 10*Y0**5*Y1/7 + 10*Y0**3*Y1*Y2**2/7 + 5*Y0*Y2**5/7 + Y1**6 + 15*Y1**5*Y2/7 + Y2**6')
A = P('X0**3*Y0**3')
F0 = add(f, g)
print('=== (a) t=0 instance, d=6, J family')
sysR = ReducedE2(F0, A, 6, 0, verbose=True)
famJ = j_family(F0, A, 6, 0, sysR.Dm, sysR.Ncols, exclude_F=True)
print('J family K =', famJ.K, 'removed F dim', getattr(famJ, 'removed_F_dim', None))
print('E1 check (rel residuals):', check_family_e1(F0, A, famJ))
Dfg = sysR.Dm.vec(sub(f, g))
a = famJ.Dbasis.T @ Dfg
print('projection residual of f-g onto family:', np.linalg.norm(famJ.Dbasis @ a - Dfg) / np.linalg.norm(Dfg))
for kz in (False, True):
    fval, rel, info = sysR.value(a, famJ, koszul=kz)
    print(f'  koszul={kz}: f={fval:.3e} rel={rel:.3e} bnorm={info["bnorm"]:.3e}')
rng = np.random.default_rng(1); ar = rng.standard_normal(famJ.K)
for kz in (False, True):
    fval, rel, info = sysR.value(ar, famJ, koszul=kz)
    print(f'  random a koszul={kz}: rel={rel:.3e}')
# gradient check
u0 = rng.standard_normal(famJ.K); u0 /= np.linalg.norm(u0)
f0, g0 = sysR.value_and_grad(u0, famJ, koszul=True)
h = 1e-6; gnum = np.zeros(famJ.K)
for k in range(5):
    e = np.zeros(famJ.K); e[k] = h
    gnum[k] = (sysR.value_and_grad(u0 + e, famJ)[0] - sysR.value_and_grad(u0 - e, famJ)[0]) / (2*h)
print('grad check (first 5):', np.c_[g0[:5], gnum[:5]])
t0 = time.time()
cands = sysR.search(famJ, n_restarts=20, seed=3, koszul=True, verbose=False)
print('search t=0 (20 restarts): best rels', [f"{c['rel']:.2e}" for c in cands[:6]], f'{time.time()-t0:.1f}s')
best = cands[0]['a']; Dv = famJ.D_vec(best)
# compare with span{f,g}: project Dv to span{f,g}
fv, gv = sysR.Dm.vec(f), sysR.Dm.vec(g); Bfg = np.linalg.qr(np.c_[fv, gv])[0]
print('best D distance from span{f,g}:', np.linalg.norm(Dv - Bfg @ (Bfg.T @ Dv)))
print('=== (b) accepted witness family-shear (d=6)')
Fs, As, qs, rec = core.load_instance('family-shear')
Wd = json.load(open('/home/kejace/devel/sextics/data/family-shear/witness.json'))
Ds = P(Wd['D'])
sysS = ReducedE2(Fs, As, 6, 0, verbose=True)
famS = full_family(Fs, As, 6, 0, verbose=True)
print('family-shear full family K =', famS.K, famS.info)
print('E1 check:', check_family_e1(Fs, As, famS))
Dsv = sysS.Dm.vec(Ds); aS = famS.Dbasis.T @ Dsv
print('projection residual of witness D onto family:', np.linalg.norm(famS.Dbasis @ aS - Dsv) / np.linalg.norm(Dsv))
for kz in (False, True):
    fval, rel, info = sysS.value(aS, famS, koszul=kz)
    print(f'  witness D koszul={kz}: rel={rel:.3e}')
for kz in (False, True):
    fval, rel, info = sysS.value(rng.standard_normal(famS.K), famS, koszul=kz)
    print(f'  random D koszul={kz}: rel={rel:.3e}')
print('e2_modp_check witness:', e2_modp_check(Fs, As, Ds))
Dr = exact_D_from_frac(sysS.Dm, [Fraction(int(x)) for x in np.round(3*rng.standard_normal(sysS.Dm.n))])
print('e2_modp_check random:', e2_modp_check(Fs, As, Dr))
t0 = time.time()
cands = sysS.search(famS, n_restarts=10, seed=5, koszul=True)
print('search family-shear (10 restarts): best rels', [f"{c['rel']:.2e}" for c in cands[:5]], f'{time.time()-t0:.1f}s')
if cands[0]['rel'] < 1e-10:
    frs = rationalize_vec(famS.D_vec(cands[0]['a']), max_den=200)
    print('rationalized:', [core.to_str(exact_D_from_frac(sysS.Dm, fr))[:120] for fr in frs])
    for fr in frs[:1]:
        print('exact check:', {k: v for k, v in exact_check_D(Fs, As, 1, exact_D_from_frac(sysS.Dm, fr)).items() if k not in ('N',)})
