"""w2s_e1spaces.py -- E1 bookkeeping (numerical SVD, both tau-parities, d = 4..12) + exact mod-p cross-checks + kernel structure.
Outputs: w2s_cache/*.npz (J complements, K_full bases, lifts), w2s_e1spaces.json, log events."""
import sys, os, time, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2s_common import *
from w2s_nvp import family

def hilbert():
    """exact dims of S_k^+/- and R_k^+/- for a smooth sextic in 6 variables (parity-graded)."""
    def pm(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): out[i + j] += x * y
        return out
    H = [1]
    for _ in range(6): H = pm(H, [1, 1, 1, 1, 1])
    Th = [1]
    for _ in range(3): Th = pm(Th, [1, 0, 1, 0, 1, 0, 1, 0, 1])
    H += [0] * 40; Th += [0] * 40
    return {k: ((H[k] + Th[k]) // 2, (H[k] - Th[k]) // 2) for k in range(30)}

def main():
    F, A, q = load()
    Fi = [vec(pderiv(F, i), 5) for i in range(NV)]
    R = hilbert()
    table = []
    t00 = time.time()
    for d in range(4, 13):
        for pD in (0, 1):
            t0 = time.time()
            Kb, rankA, cols = K_full(F, A, Fi, d, pD)
            Jb, rJ = J_basis(Fi, d, pD)
            # numerical check J subset K_full
            resJ = np.linalg.norm(Jb - Kb @ (Kb.T @ Jb)) / max(np.linalg.norm(Jb), 1e-300) if rJ else 0.0
            rec = dict(d=d, pD=pD, dimS=int(len(cols)), dimJ_num=int(rJ), dimJ_exact=int(len(cols) - R[d][pD]),
                       dimR_exact=int(R[d][pD]), dimKfull_num=int(Kb.shape[1]), kernel_modJ_num=int(Kb.shape[1] - rJ),
                       rankA_num=int(rankA), dimR_target_exact=int(R[d + 6][(pD + 1) % 2]),
                       min_kernel_by_count=int(max(0, R[d][pD] - R[d + 6][(pD + 1) % 2])),
                       J_in_Kfull_resid=float(resJ), seconds=round(time.time() - t0, 1))
            table.append(rec)
            log_event(event='e1-space', **rec)
            print(json.dumps(rec), flush=True)
            with open(os.path.join(LANE, 'w2s_e1spaces.json'), 'w') as f:
                json.dump(table, f, indent=1)
    print(f'numerical table done in {time.time()-t00:.0f}s', flush=True)
    # families with lifts for d = 10 (both parities), kernel-only and full
    for pD in (0, 1):
        for kind in ('ker', 'kfull'):
            t0 = time.time()
            Dmat, N0, P, lab = family(F, A, 10, pD, kind)
            log_event(event='family-built', d=10, pD=pD, kind=kind, dim=int(Dmat.shape[0]), label=lab, seconds=round(time.time() - t0, 1))
            print('family', lab, time.time() - t0, flush=True)
    # kernel structure: coordinate subspaces
    E_cache = {}
    def subspaces(d):
        E = E_of(d)
        X = E[:, :3].sum(1); Y = E[:, 3:].sum(1)
        out = {'X-only': Y == 0, 'Y-only': X == 0, 'vars X0X1X2Y0': (E[:, 4] + E[:, 5]) == 0,
               'vars X0Y0Y1Y2': (E[:, 1] + E[:, 2]) == 0, 'vars X0Y0': (E[:, 1] + E[:, 2] + E[:, 4] + E[:, 5]) == 0,
               'vars X0X1Y0Y1': (E[:, 2] + E[:, 5]) == 0, 'div X0': E[:, 0] >= 1, 'div Y0': E[:, 3] >= 1,
               'div X0Y0': (E[:, 0] >= 1) & (E[:, 3] >= 1), 'div X0^2': E[:, 0] >= 2, 'div Y0^2': E[:, 3] >= 2,
               'div X0^2Y0^2': (E[:, 0] >= 2) & (E[:, 3] >= 2), 'div X0^3Y0^3': (E[:, 0] >= 3) & (E[:, 3] >= 3)}
        for k in range(2, 7):
            out[f'(X)^{k}'] = X >= k; out[f'(Y)^{k}'] = Y >= k
            out[f'(X0,Y0)^{k}'] = (E[:, 0] + E[:, 3]) >= k
            out[f'(X1,X2,Y1,Y2)^{k}'] = (E[:, 1] + E[:, 2] + E[:, 4] + E[:, 5]) >= k
            out[f'(X1,X2)^{k}'] = (E[:, 1] + E[:, 2]) >= k
            out[f'(Y1,Y2)^{k}'] = (E[:, 4] + E[:, 5]) >= k
            out[f'(X0,X1,Y0,Y1)^{k}'] = (E[:, 0] + E[:, 1] + E[:, 3] + E[:, 4]) >= k
        return out
    struct = []
    for d in (10, 11, 12):
        for pD in (0, 1):
            Kb, rankA, cols = K_full(F, A, Fi, d, pD)
            Jb, rJ = J_basis(Fi, d, pD)
            dimK = Kb.shape[1]
            for name, mask in subspaces(d).items():
                notV = np.nonzero(~mask)[0]
                def dim_inter(B):
                    if B.shape[1] == 0: return 0
                    M = B[notV, :]
                    if M.shape[0] == 0: return B.shape[1]
                    s = np.linalg.svd(M, compute_uv=False)
                    r = int(np.sum(s > 1e-9 * max(s[0], 1e-300)))
                    return B.shape[1] - r
                kV = dim_inter(Kb); jV = dim_inter(Jb)
                dimV = int(mask.sum() if pD is None else np.sum(mask & (par_of(d) == pD)))
                rec = dict(d=d, pD=pD, subspace=name, dimV=dimV, dim_Kfull_cap_V=int(kV), dim_J_cap_V=int(jV), new_nonjacobian=int(kV - jV))
                struct.append(rec)
                if kV - jV:
                    print(json.dumps(rec), flush=True)
            log_event(event='kernel-structure', d=d, pD=pD, dimKfull=int(dimK), dimJ=int(rJ),
                      nonzero=[(r['subspace'], r['dimV'], r['new_nonjacobian']) for r in struct if r['d'] == d and r['pD'] == pD and r['new_nonjacobian']])
            with open(os.path.join(LANE, 'w2s_kernel_structure.json'), 'w') as f:
                json.dump(struct, f, indent=1)
    print('structure done', flush=True)
    # exact mod-p cross-checks with the engine (single-threaded FLINT)
    for d in range(4, 11):
        for pD in (0, 1):
            t0 = time.time()
            try:
                out = jacobian.mult_kernel(F, A, d, parity=pD)
                rec = dict(event='e1-exact-modp', d=d, pD=pD, dim_kernel_modJ=int(out['dim']), kernel_dim=int(out['kernel_dim']),
                           jacobian_dim=int(out['jacobian_dim']), ring_dim=int(out['ring_dim']), p=int(out['p']), seconds=round(time.time() - t0, 1))
                if d == 10:
                    with open(os.path.join(LANE, f'w2s_exact_kernel_d10_par{pD}.json'), 'w') as f:
                        json.dump({'p': out['p'], 'basis': [{' '.join(map(str, k)): v for k, v in b.items()} for b in out['basis']]}, f)
            except Exception as ex:
                rec = dict(event='e1-exact-modp', d=d, pD=pD, error=str(ex)[:200], seconds=round(time.time() - t0, 1))
            log_event(**rec)
            print(json.dumps(rec), flush=True)

if __name__ == '__main__':
    main()
