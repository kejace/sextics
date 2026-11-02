"""Wave-2 job: float basis of the E1-admissible denominator space (kernel + J_d) with particular E1 solutions.

E1: A D - sum_i N_i F_i = 0.  Columns [A m (m in S_d^{pD}) | -m F_i (N-monomials)]; the null space projected
to the D-coordinates is the admissible space (dimension known exactly mod p: 1044 even / 1008 odd at d = 10).
Output npz: Dfull (nD x K), Nfull (nNtot x K) with N_k(F) = A D_k, singular-value gap diagnostics.

usage: python3 w2_k10prep.py --d 10 --par 0 [--threads 8] [--K 1044]
"""
import sys, os, json, time, argparse
ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, default=10)
ap.add_argument('--par', type=int, default=0)
ap.add_argument('--threads', type=int, default=8)
ap.add_argument('--K', type=int, default=None, help='expected admissible dimension (else from gap)')
ap.add_argument('--tag', default='k10prep')
args = ap.parse_args()
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = str(args.threads)
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *

d, pD = args.d, args.par
t0 = time.time()
pars = [(pA - pD - (0 if i < 3 else 1)) % 2 for i in range(NV)]
Ds = Space(d, pD)
Ns = [Space(d + 1, pars[i]) for i in range(NV)]
Ts = Space(d + 6, 1)
cols = [mul_matrix(A, Ds, Ts)] + [-mul_matrix(Fi[i], Ns[i], Ts) for i in range(NV)]
M = sp.hstack(cols).tocsr()
nD = Ds.n; nN = sum(N.n for N in Ns)
print(f'[k10prep] d={d} par={pD} E1 matrix {M.shape} nnz={M.nnz} nD={nD} nN={nN} ({time.time()-t0:.1f}s)', file=sys.stderr, flush=True)
Md = M.toarray()
t1 = time.time()
U, s, Vt = sla.svd(Md, full_matrices=True, check_finite=False, lapack_driver='gesdd')
t2 = time.time()
tol = s[0] * 1e-11
rank = int(np.sum(s > tol))
null = Vt[rank:].T            # (nD + nN) x nz
nz = null.shape[1]
print(f'[k10prep] SVD {t2-t1:.1f}s rank={rank} nullity={nz} s[rank-1]={s[rank-1]:.3e} s[rank]={(s[rank] if rank < len(s) else 0):.3e}', file=sys.stderr, flush=True)
Zd = null[:nD]; Zn = null[nD:]
Ud, sd, Vdt = sla.svd(Zd, full_matrices=False, check_finite=False)
K = args.K
if K is None:
    gaps = sd[:-1] / np.maximum(sd[1:], 1e-300)
    K = int(np.argmax(gaps)) + 1
Dfull = Ud[:, :K]
Nfull = Zn @ (Vdt[:K].T / sd[:K])
# verify E1 for the basis
chk = Md @ np.vstack([Dfull, Nfull])
e1err = float(np.abs(chk).max())
# verify that the J-part is inside: J_d generators m F_i should be in the span of Dfull
Wss = [Space(d - 5, (pD - (0 if i < 3 else 1)) % 2) for i in range(NV)]
Jcols = sp.hstack([mul_matrix(Fi[i], Wss[i], Ds) for i in range(NV)]).toarray()
proj = Dfull @ (Dfull.T @ Jcols)
jerr = float(np.linalg.norm(proj - Jcols) / np.linalg.norm(Jcols))
out = f'{WD}/w2_kernelfam_d{d}_p{pD}.npz'
np.savez(out, Dfull=Dfull, Nfull=Nfull, sd=sd[:K + 5], s_e1=s[max(0, rank - 3):rank + 3])
info = {'event': 'k10prep', 'd': d, 'par': pD, 'nD': nD, 'nN': nN, 'rows': Ts.n, 'rank_E1': rank, 'nullity': nz,
        'K': K, 'sd_K': float(sd[K - 1]), 'sd_K+1': float(sd[K]) if K < len(sd) else None,
        'gap_ratio': float(sd[K - 1] / sd[K]) if K < len(sd) else None, 'e1_max_err': e1err, 'J_in_span_rel_err': jerr,
        'file': out, 'seconds': round(time.time() - t0, 1)}
jlog(f'{WD}/{args.tag}.jsonl', info)
print(json.dumps(info), flush=True)
