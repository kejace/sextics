"""Self-test of the normalized-E2 float machinery against accepted witnesses and known identities."""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
import numpy as np

def lstsq_reduced(prob, Dv, Nv, cond=1e-10):
    B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
    Br = prob.reduce(B); Mr = prob.reduce(M)
    c, *_ = sla.lstsq(Mr, -Br, cond=cond, lapack_driver='gelsd', check_finite=False)
    r = Br + Mr @ c
    return float(np.linalg.norm(r) / max(np.linalg.norm(Br), 1e-300)), float(np.linalg.norm(Br) / max(np.linalg.norm(B), 1e-300)), c

def test_witness(iid):
    print(f'=== witness test {iid}', flush=True)
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    print('  verify_witness', ok, info)
    N2 = normalize_lift(F, A, D, N)
    pD = parity_of(D)
    prob = E2Float(F, A, d, pD)
    Dv = prob.SD.vec(D); Nv = prob.Nvec(N2)
    B = prob.Bvec(Dv, Nv)
    Bex = prob.Srow.vec(B_of(D, N2))
    print('  |B float - B exact| / |B| =', np.linalg.norm(B - Bex) / np.linalg.norm(Bex))
    print('  |rho B| / |B| (should be ~0: B in F*S) =', np.linalg.norm(prob.reduce(B)) / np.linalg.norm(B))
    rng = np.random.default_rng(1)
    c = rng.standard_normal(prob.n_c)
    Nk = Nv + prob.Kmap @ c
    rel, _, cc = lstsq_reduced(prob, Dv, Nk)
    print('  perturbed by Koszul field: rel after varpro =', rel, ' |c_rec + c| =', np.linalg.norm(cc + c) / np.linalg.norm(c))
    Nr = Nv + 1e-3 * rng.standard_normal(prob.nN) * np.linalg.norm(Nv) / np.sqrt(prob.nN)
    rel2, _, _ = lstsq_reduced(prob, Dv, Nr)
    print('  perturbed by random non-Koszul field: rel =', rel2)
    fam = E1Family(prob)
    co = fam.Dbasis @ Dv
    print('  witness D in family: fit =', np.linalg.norm(fam.Dbasis.T @ co - Dv) / np.linalg.norm(Dv), ' K =', fam.K, ' J_d dim =', fam.J_part_dim())
    N0 = fam.Nbasis.T @ co
    rel3, frac3, _ = lstsq_reduced(prob, Dv, N0)
    print('  varpro residual at witness point with family lift: rel =', rel3)
    # random family point should be far from zero
    a = rng.standard_normal(fam.K)
    rel4, _, _ = lstsq_reduced(prob, fam.Dbasis.T @ a, fam.Nbasis.T @ a)
    print('  random family point: rel =', rel4)
    # gradient check of f(a) = |r|^2 via adjoint
    def fval(a):
        Dv = fam.Dbasis.T @ a; Nv = fam.Nbasis.T @ a
        B = prob.Bvec(Dv, Nv); M = prob.koszul_matrix(Dv)
        Br = prob.reduce(B); Mr = prob.reduce(M)
        c, *_ = sla.lstsq(Mr, -Br, cond=1e-10, lapack_driver='gelsd', check_finite=False)
        r = Br + Mr @ c
        w = prob.reduce_T(r)
        Ntot = Nv + prob.Kmap @ c
        gD, gN = prob.adjoint_B(w, Dv, Ntot)
        g = 2.0 * (fam.Dbasis @ gD + fam.Nbasis @ gN)
        return float(r @ r), g
    f0, g0 = fval(a)
    h = rng.standard_normal(fam.K); h /= np.linalg.norm(h)
    eps = 1e-6
    fp, _ = fval(a + eps * h); fm, _ = fval(a - eps * h)
    print('  gradient check: analytic', g0 @ h, ' numeric', (fp - fm) / (2 * eps))
    return prob, fam

def test_lane():
    print('=== lane identities', flush=True)
    F, A, q = load()
    st = structure(F)
    RX = euler_part(F, [0, 1, 2])
    X = [{tuple(1 if j == i else 0 for j in range(6)): Fraction(1)} for i in range(3)]
    N0 = [mul(A, X[0]), mul(A, X[1]), mul(A, X[2]), {}, {}, {}]
    B = B_of(RX, N0)
    target = scale(mul(A, st['P']), 8 * st['t'])
    print('  B_{eps_X F}(A X) == 8 t A L^2 M :', sub(B, target) == {})
    prob = E2Float(F, A, 6, 0)
    Dv = prob.SD.vec(RX); Nv = prob.Nvec(N0)
    Bf = prob.Bvec(Dv, Nv)
    print('  float B matches exact:', np.linalg.norm(Bf - prob.Srow.vec(B)) / np.linalg.norm(Bf))
    rel, frac, _ = lstsq_reduced(prob, Dv, Nv)
    print(f'  full normalized E2 for D = eps_X F: rel residual (reduced) = {rel:.3e}; |rho B|/|B| = {frac:.3e}')
    relf, _, _ = prob.residual_full(Dv, Nv)
    print(f'  un-reduced residual (with C block): {relf:.3e}')
    fam = E1Family(prob)
    print('  d=6 family K =', fam.K, '(expect 18 = dim J_6^+), J part', fam.J_part_dim(), ' spurious dims', fam.spurious_subspace().shape[1])
    co = fam.Dbasis @ Dv
    print('  eps_X F in family fit', np.linalg.norm(fam.Dbasis.T @ co - Dv) / np.linalg.norm(Dv))
    # generic rank of the reduced Koszul matrix
    Mr = prob.reduce(prob.koszul_matrix(Dv))
    s = sla.svdvals(Mr)
    rk = int((s > 1e-9 * s[0]).sum())
    print(f'  reduced Koszul matrix {Mr.shape}: numerical rank {rk}, nullity {Mr.shape[1]-rk}, smallest kept {s[rk-1]:.2e}, next {s[rk] if rk < len(s) else 0:.2e}')

if __name__ == '__main__':
    which = sys.argv[1:] or ['coupled-pair-rational-d4-strict', 'lane']
    for w in which:
        t0 = time.time()
        if w == 'lane':
            test_lane()
        else:
            test_witness(w)
        print(f'  ({time.time()-t0:.1f}s)', flush=True)
