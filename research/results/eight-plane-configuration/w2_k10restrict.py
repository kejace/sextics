"""Restrict the float admissible family (w2_kernelfam_d10_p{par}.npz: Dfull, Nfull) to a structured subspace:
denominators vanishing to order k along all eight planes, i.e. D in (Q0,Q1,Q2)^k * S_{10-2k}.
Computes the intersection in floats (SVD of [Dfull | -ideal basis]) and writes w2_kernelfam_d10_p{par}_Q{k}.npz
with Dfull_r = Dfull @ X (nD x K_r) and Nfull_r = Nfull @ X (so the particular E1 solutions restrict consistently).

usage: python3 w2_k10restrict.py --par 0 --k 2 [--threads 8]
"""
import sys, os, json, time, argparse, itertools
ap = argparse.ArgumentParser()
ap.add_argument('--par', type=int, default=0)
ap.add_argument('--k', type=int, default=2)
ap.add_argument('--d', type=int, default=10)
ap.add_argument('--threads', type=int, default=8)
args = ap.parse_args()
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = str(args.threads)
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *

d, par, k = args.d, args.par, args.k
z = np.load(f'{WD}/w2_kernelfam_d{d}_p{par}.npz')
Dfull, Nfull = z['Dfull'], z['Nfull']
Ds = Space(d, par)
Qs = [core.parse('X0**2+Y0**2'), core.parse('X1**2+Y1**2'), core.parse('X2**2+Y2**2')]
# ideal basis: products of k Q's (with repetition) times monomials of degree d-2k, parity par
cols = []
for combo in itertools.combinations_with_replacement(range(3), k):
    P = {(0,) * 6: Fraction(1)}
    for j in combo:
        P = mul(P, Qs[j])
    src = Space(d - 2 * k, par)
    cols.append(mul_matrix(P, src, Ds))
Ideal = sp.hstack(cols).toarray()
# orthonormal basis of the ideal subspace
Ui, si, _ = np.linalg.svd(Ideal, full_matrices=False)
ri = int(np.sum(si > 1e-10 * si[0]))
Ui = Ui[:, :ri]
# intersection: a with Dfull a in span(Ui): nullspace of (I - Ui Ui^T) Dfull
Mres = Dfull - Ui @ (Ui.T @ Dfull)
X = sla.null_space(Mres, rcond=1e-9)
Kr = X.shape[1]
Dr = Dfull @ X; Nr = Nfull @ X
# orthonormalize the restricted D-basis (keep N consistent)
Ud, sd, Vdt = np.linalg.svd(Dr, full_matrices=False)
rk = int(np.sum(sd > 1e-10 * sd[0]))
T = Vdt[:rk].T / sd[:rk]
Dr2 = Dr @ T; Nr2 = Nr @ T
out = f'{WD}/w2_kernelfam_d{d}_p{par}_Q{k}.npz'
np.savez(out, Dfull=Dr2, Nfull=Nr2)
info = {'event': 'k10restrict', 'd': d, 'par': par, 'k': k, 'ideal_dim': ri, 'K_full': Dfull.shape[1], 'K_restricted': rk,
        'residual_of_restriction': float(np.linalg.norm(Dr2 - Ui @ (Ui.T @ Dr2)) / max(np.linalg.norm(Dr2), 1e-300)), 'file': out}
jlog(f'{WD}/k10prep.jsonl', info)
print(json.dumps(info), flush=True)
