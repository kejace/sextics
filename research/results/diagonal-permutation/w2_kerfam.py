"""Float basis of the E1-admissible denominators at degree d (Y-even) with E1 lifts, by SVD of the E1 matrix
[ A*m (m in S_d^+) | -F_i*m' (N-monomials) ] over the reals.  Output: w2-kerfam-d<d>.npz with
  Dmons (n x 6), Dvecs (K x n) float coefficient vectors spanning a complement of J_d^+ inside the admissible space,
  Nmons_i, Nvecs (K x nN) lifts with N(F) = A D (to float accuracy), gap info.
usage: python3 w2_kerfam.py d [pD]
"""
import sys, os, time, json
import numpy as np, scipy.linalg as sla
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from hodge.core import monomials, NV, pderiv, shift, mul
LANE = os.path.dirname(os.path.abspath(__file__))
d = int(sys.argv[1]); pD = int(sys.argv[2]) if len(sys.argv) > 2 else 0
F, A, q, rec = core.load_instance('diagonal-permutation')
APRIME = sys.argv[3] if len(sys.argv) > 3 else ''
if APRIME:
    A = core.parse(json.load(open(os.path.join(LANE, APRIME)))['A'])
    print(f'[kerfam] target overridden by {APRIME}', flush=True)
t0 = time.time()
PAR = [0, 0, 0, 1, 1, 1]
Dmons = monomials(d, pD); nD = len(Dmons)
pars = fixed_d.numerator_parities({Dmons[0]: 1}, A)
Ncols = [(i, m) for i in range(NV) for m in monomials(d + 1, pars[i])]
nN = len(Ncols)
rows = {e: j for j, e in enumerate(monomials(d + 6, (1 + pD) % 2))}
R = len(rows); C = nD + nN
M = np.zeros((R, C))
for j, m in enumerate(Dmons):
    for e, v in shift(A, m).items():
        M[rows[e], j] = float(v)
Fi = [pderiv(F, i) for i in range(NV)]
for j, (i, m) in enumerate(Ncols):
    for e, v in shift(Fi[i], m).items():
        M[rows[e], nD + j] = -float(v)
print(f'[d={d}] E1 matrix {R} x {C} built in {time.time()-t0:.1f}s', flush=True)
U, s, Vt = sla.svd(M, full_matrices=True, lapack_driver='gesdd', check_finite=False)
t1 = time.time()
tol = 1e-9 * s[0]
rank = int(np.sum(s > tol))
gap = (s[rank - 1], s[rank] if rank < len(s) else 0.0)
null = Vt[rank:].T                       # C x nullity
nullity = null.shape[1]
print(f'[d={d}] svd {t1-t0:.1f}s rank={rank} nullity={nullity} gap=({gap[0]:.3e}, {gap[1]:.3e})', flush=True)
# admissible D-space = projection of the null space to D-coordinates
Na = null[:nD, :]
Ua, sa, Vat = sla.svd(Na, full_matrices=False, check_finite=False)
tola = 1e-9 * sa[0]
ra = int(np.sum(sa > tola))
print(f'[d={d}] admissible dim = {ra} (singular gap {sa[ra-1]:.3e} / {sa[ra] if ra < len(sa) else 0:.3e})', flush=True)
# basis of the admissible space: columns of Ua[:, :ra]; lifts: null @ Vat[k,:]/sa[k]
Vadm = Ua[:, :ra]                                        # nD x ra
Lift = null @ (Vat[:ra, :].T / sa[:ra])                  # C x ra  (D-part = Vadm)
assert np.allclose(Lift[:nD], Vadm, atol=1e-8)
# Jacobian part J_d^+ (exact generators m F_i) as float vectors
Didx = {m: j for j, m in enumerate(Dmons)}
Jgens = []
for i in range(NV):
    for m in monomials(d - 5, (pD + PAR[i]) % 2):
        v = np.zeros(nD)
        for e, c in shift(Fi[i], m).items():
            v[Didx[e]] = float(c)
        Jgens.append(v)
Jm = np.array(Jgens).T                                   # nD x 756
Qj, _ = np.linalg.qr(Jm)
Cperp = Vadm - Qj @ (Qj.T @ Vadm)
Uc, sc, Vct = sla.svd(Cperp, full_matrices=False, check_finite=False)
rk = int(np.sum(sc > 1e-9 * max(1.0, sc[0])))
print(f'[d={d}] complement of J in admissible: dim {rk} (J rank {np.linalg.matrix_rank(Jm)}); gap {sc[rk-1]:.3e}/{sc[rk] if rk < len(sc) else 0:.3e}', flush=True)
alpha = Vct[:rk, :].T                                    # ra x rk : c_k = Vadm @ alpha_k
Dk = (Vadm @ alpha).T                                    # rk x nD
Nk = (Lift[nD:, :] @ alpha).T                            # rk x nN
# verify E1 residual for a few
res = []
for k in range(min(3, rk)):
    D = {Dmons[j]: Dk[k, j] for j in range(nD) if abs(Dk[k, j]) > 1e-15}
    N = [{} for _ in range(NV)]
    for t, (i, m) in enumerate(Ncols):
        if abs(Nk[k, t]) > 1e-15: N[i][m] = Nk[k, t]
    lhs = {}
    for i in range(NV):
        for e, v in mul(N[i], Fi[i]).items(): lhs[e] = lhs.get(e, 0) + v
    rhs = mul(A, D)
    err = max(abs(lhs.get(e, 0) - rhs.get(e, 0)) for e in set(lhs) | set(rhs))
    res.append(err)
print(f'[d={d}] E1 residuals of first members: {res}', flush=True)
np.savez_compressed(os.path.join(LANE, f'w2-kerfam-d{d}-p{pD}' + (('-' + APRIME.replace('.json','')) if APRIME else '') + '.npz'), Dmons=np.array(Dmons), Dvecs=Dk, Ncols_i=np.array([c[0] for c in Ncols]),
                    Ncols_m=np.array([c[1] for c in Ncols]), Nvecs=Nk, sv_gap=np.array(gap), adm_dim=ra, ker_dim=rk)
rec = {'exp': 'kerfam-float', 'd': d, 'pD': pD, 'Aprime': APRIME, 'rows': R, 'cols': C, 'rank': rank, 'nullity': nullity, 'gap': [float(g) for g in gap],
       'admissible_dim': int(ra), 'kernel_mod_J_dim': int(rk), 'e1_res': [float(x) for x in res], 'seconds': round(time.time() - t0, 1)}
with open(os.path.join(LANE, 'w2-kerfam.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
print(json.dumps(rec), flush=True)
