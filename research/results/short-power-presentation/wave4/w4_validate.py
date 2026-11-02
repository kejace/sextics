import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import fam_vert, restricted_kernel
F, A, q, rec = core.load_instance(IID)
p = core.PRIMES[0]
Fp = reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
def comb_(P, Q, c):  # P - c Q
    return core.sub(P, core.scale(Q, c, p), p)
chk = [reduce_T(comb_(Fi[1], Fi[0], 2), p), reduce_T(comb_(Fi[2], Fi[0], 3), p),
       reduce_T(comb_(Fi[4], Fi[3], 3), p), reduce_T(comb_(Fi[5], Fi[3], 5), p)]
print('I-prime relations reduce to zero in T:', [len(x) for x in chk])
print('F0 terms', len(Fi[0]), 'F3 terms', len(Fi[3]))
# Hilbert: engine ranks rank J_12^- = 2214 (target odd), J_12^+ = 2223
print('dimS(12,1)-dimR(12,1) =', dimS(12,1)-dimR(12,1), '(engine 2214)   dimS(12,0)-dimR(12,0) =', dimS(12,0)-dimR(12,0), '(engine 2223)')
print('dimR total by degree 0..24:', [dimR(n,0)+dimR(n,1) for n in range(25)])
print('tbasis counts vs dimT:', [(n, len(tbasis(n,0)), dimT(n,0), len(tbasis(n,1)), dimT(n,1)) for n in (5, 12, 13, 18)])
# unrestricted E1 kernels d=5,6 (engine: kernel 3 / 18, kernel mod J = 0)
for d in (5, 6):
    for par in (0, 1):
        fam = [(str(m), {m: 1}) for m in monomials(d, par)]
        t0 = time.time(); vecs, info = e1_kernel(F, A, fam, d, par, p)
        vj, infoj = in_J_family(F, fam, d, par, p)
        print('unrestricted d=%d par=%d: fam %d E1-kernel %d  (J_d cap fam: %d)  hilbert_ok=%s ring_dim=%d/%d  %.1fs' % (
            d, par, len(fam), info['kernel_dim'], infoj['kernel_dim'], info['hilbert_ok'], info['ring_dim'], info['ring_dim_pred'], time.time()-t0))
        log({'what': 'validate-unrestricted-E1', 'd': d, 'parD': par, 'prime': p, 'kernel_dim': info['kernel_dim'], 'J_d_dim': infoj['kernel_dim'],
             'hilbert_ok': info['hilbert_ok'], 'engine_reference': {5: 3, 6: 18}[d], 'file': 'w4_validate.py'})
# vertical family d=6: T-method vs wave-2 restricted_kernel (big nullspace)
for par in (0, 1):
    fam = fam_vert(6, par)
    vecs, info = e1_kernel(F, A, fam, 6, par, p)
    t0 = time.time(); v2, info2 = restricted_kernel(F, A, fam, 6, p, verbose=False); t1 = time.time() - t0
    print('vert d=6 par=%d: T-method kernel %d, wave-2 nullspace kernel %d (%.0fs); T vecs %s' % (par, info['kernel_dim'], info2['kernel_dim'], t1, [sym(v, p) for v in vecs]))
    log({'what': 'validate-vert-d6-vs-restricted_kernel', 'd': 6, 'parD': par, 'prime': p, 'T_kernel': info['kernel_dim'], 'w2_kernel': info2['kernel_dim'], 'file': 'w4_validate.py'})
