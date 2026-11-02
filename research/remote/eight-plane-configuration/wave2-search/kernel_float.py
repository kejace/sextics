"""Float E1 kernel at degree d, parity parD: orthonormal D-basis X (nD x k) of {D : A D in J(F)_{d+6}} and exact-in-float lifts
Y with M_D X + M_N Y = 0 (N = lifts, sum_i N_i F_i = A D).  Saves wave2-search/kernel_float_d{d}_par{parD}.npz with
D (k x nD), N0..N5 (k x nN_i), residual norms.  Method: SVD of M_N (left null space Z), kernel = null(Z^T M_D), lifts by pseudo-inverse."""
import sys, time, json, numpy as np, scipy.linalg as sla
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-search')
from w2ep import *
d, parD = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
F, A, q = load()
parA = core.parity_of(A)
SD = space(d, parD); pars = numerator_parities(parD, parA); SN = [space(d + 1, pars[i]) for i in range(NV)]
SR = space(d + 6, (parA + parD) % 2)
Ak, Av = poly_arrays(A)
MD = mult_matrix(Ak, Av, SD, SR).toarray()
Fi = [pderiv(F, i) for i in range(NV)]
blocks = []
for i in range(NV):
    ks, vs = poly_arrays(Fi[i]); blocks.append(mult_matrix(ks, vs, SN[i], SR).toarray())
MN = np.hstack(blocks); del blocks
print(f'[kernel_float] d={d} parD={parD} rows={SR.n} nD={SD.n} nN={MN.shape[1]} built {time.time()-t0:.1f}s', flush=True)
t1 = time.time()
U, s, Vt = sla.svd(MN, full_matrices=False, lapack_driver='gesdd', check_finite=False, overwrite_a=False)
r = int(np.sum(s > 1e-10 * s[0]))
print(f'[kernel_float] svd(M_N) {time.time()-t1:.1f}s rank {r} (codim {SR.n - r}); s[0]={s[0]:.3e} s[r-1]={s[r-1]:.3e} s[r]={s[r] if r < len(s) else 0:.3e}', flush=True)
Z = U[:, r:]                       # left null space of M_N (orthonormal), dim = codim of J_{d+6}
C = Z.T @ MD                       # codim x nD
X = sla.null_space(C, rcond=1e-10)  # nD x k
k = X.shape[1]
# lifts: Y = - V_r S_r^-1 U_r^T M_D X
Y = -(Vt[:r, :].T @ ((U[:, :r].T @ (MD @ X)) / s[:r, None]))
res = np.linalg.norm(MD @ X + MN @ Y, axis=0)
print(f'[kernel_float] kernel dim {k}; max E1 residual {res.max():.3e} (|M_D X| ~ {np.linalg.norm(MD @ X, axis=0).max():.3e}); total {time.time()-t0:.1f}s', flush=True)
out = {'D': X.T.copy()}
off = 0
for i in range(NV):
    out[f'N{i}'] = Y[off:off + SN[i].n, :].T.copy(); off += SN[i].n
out['e1_residual'] = res; out['svals_MN'] = s
np.savez_compressed(f'{SD_DIR}/kernel_float_d{d}_par{parD}.npz', **out)
log_batch({'batch': 'kernel-float', 'd': d, 'parity': parD, 'rows': SR.n, 'nD': SD.n, 'nN': int(MN.shape[1]), 'rank_MN': r, 'codim_J': SR.n - r,
           'kerdim': k, 'max_e1_residual': float(res.max()), 'seconds': round(time.time() - t0, 1)})
