"""Hessian relaxation at the Jacobian-type degree d=5 (exact, mod two primes) and the J-type relaxation at d=6.

d=5: every E1-admissible D of degree 5 is D = sum_i w_i F_i (constant w).  With the canonical lift N0 = (A/2) w,
B_D(N0) = (1/2)[D div(A w) - A w.grad D] and modulo J(F) only the Hessian term survives:
    [B_D(N0)] = -(1/2) sum_{i,j} w_i w_j [A F_ij]  in R_16 = (S/J)_16.
E2 forces [B_D(N0)] = 0 (im T_D and F^2 Q lie in J).  If the 21 classes [A F_ij] (i<=j) are linearly independent in R_16
then w_i w_j = 0 for all i,j, hence w = 0: no degree-5 denominator exists at all (any parity, any symmetry).
Independence mod p implies independence over Q, so the conclusion is rigorous in characteristic 0.
Classes of different G-character are automatically independent, so the test is done per character sector of R_16.

d=6 (J-part): D = sum_i l_i F_i with l_i linear forms (36 parameters, mixed sectors).  Relaxation: sum_{i,j} [A F_ij l_i l_j] = 0
in R_18.  We compute the rank of the linear map (S_2)^{21} -> R_18, c_{ij,m} -> [A F_ij m], the dimension of its kernel,
and run a cheap numerical search for w with Phi(w) = 0 on the sphere (pure polynomial objective, no inner solve)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import w2lib as W
import numpy as np
from flint import nmod_mat
from fractions import Fraction
from hodge import core

VARS = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']


def sector_nf_rank(S, deg, key, polys, p):
    """rank of the classes of `polys` (all of character key, degree deg) in (R_deg)_key mod p, plus dim (R_deg)_key"""
    mons = S.mons(deg, key)
    idx = {tuple(int(x) for x in m): j for j, m in enumerate(mons)}
    n = len(mons)
    jb = W.jac_basis(S, deg, key)
    rows = []
    for c in jb:
        row = [0] * n
        for e, v in core.reduce_mod(c, p).items():
            row[idx[e]] = v
        rows.append(row)
    rJ = 0
    if rows:
        Jm = nmod_mat(len(rows), n, [x for r in rows for x in r], p)
        rJ = Jm.rank()
    prow = []
    for P in polys:
        row = [0] * n
        for e, v in core.reduce_mod(P, p).items():
            row[idx[e]] = v
        prow.append(row)
    allm = nmod_mat(len(rows) + len(prow), n, [x for r in rows + prow for x in r], p)
    rtot = allm.rank()
    return {'dimS': n, 'dimJ': rJ, 'dimR': n - rJ, 'rank_classes': rtot - rJ, 'n_classes': len(polys)}


def d5_test(p):
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    Fi = [core.pderiv(F, i) for i in range(6)]
    Fij = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    groups = {}
    for i in range(6):
        for j in range(i, 6):
            P = core.mul(A, Fij[i][j])
            if not P:
                groups.setdefault(('zero',), []).append((i, j))
                continue
            key = S.key(P)
            groups.setdefault(key, []).append((i, j))
    out = {'p': p, 'sectors': []}
    total_rank = 0
    n_classes = 0
    zero_pairs = []
    for key, pairs in groups.items():
        if key == ('zero',):
            zero_pairs = pairs
            continue
        polys = [core.mul(A, Fij[i][j]) for (i, j) in pairs]
        r = sector_nf_rank(S, 16, key, polys, p)
        r['pairs'] = [f'{VARS[i]}{VARS[j]}' for (i, j) in pairs]
        out['sectors'].append(r)
        total_rank += r['rank_classes']; n_classes += len(pairs)
    out['total_rank'] = total_rank; out['n_classes'] = n_classes; out['zero_products'] = [f'{VARS[i]}{VARS[j]}' for (i, j) in zero_pairs]
    out['independent'] = (total_rank == n_classes) and not zero_pairs
    return out


def d6_jtype(p, n_restarts=200, seed=0):
    """relaxation for D = sum_i l_i F_i: Phi(w) = sum_{i<=j} [A F_ij l_i l_j] in R_18; rank of the linear map and numeric zeros."""
    from scipy.optimize import minimize
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    Fi = [core.pderiv(F, i) for i in range(6)]
    Fij = [[core.pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    lin = [tuple(1 if t == k else 0 for t in range(6)) for k in range(6)]
    # classes [A F_ij x_a x_b] for all i<=j, a<=b, grouped by sector of R_18; build a float normal-form basis per sector
    # (over Q we use a modular rank to count; for the numeric search we use the real quotient map via least squares on J)
    t0 = time.time()
    items = []  # (i,j,a,b, poly)
    for i in range(6):
        for j in range(i, 6):
            for a in range(6):
                for b in range(a, 6):
                    m = {tuple(lin[a][t] + lin[b][t] for t in range(6)): Fraction(1)}
                    P = core.mul(core.mul(A, Fij[i][j]), m)
                    if P:
                        items.append((i, j, a, b, P))
    bykey = {}
    for it in items:
        bykey.setdefault(S.key(it[4]), []).append(it)
    rank_info = []
    # numeric quotient projectors per sector: Pperp (R_18 complement of J in the sector)
    projs = {}
    for key, its in bykey.items():
        r = sector_nf_rank(S, 18, key, [it[4] for it in its], p)
        rank_info.append({'key': key, **r})
        mons = S.mons(18, key)
        idx = {tuple(int(x) for x in m): j for j, m in enumerate(mons)}
        Jv = W.vecs_from_dicts(mons, W.jac_basis(S, 18, key))
        Jb = W.orth_basis(Jv) if Jv.size else np.zeros((len(mons), 0))
        V = W.vecs_from_dicts(mons, [it[4] for it in its])
        Vp = V - Jb @ (Jb.T @ V) if Jb.shape[1] else V   # projection to J^perp (coordinates of the class)
        projs[key] = (its, Vp)
    tot_rank = sum(r['rank_classes'] for r in rank_info)
    dimR = sum(r['dimR'] for r in rank_info)
    # numeric search: w in R^36 (w[i][a] coefficient of x_a in l_i), Phi(w) = sum over items Vp[:, t] * w_i,a w_j,b (sym factor)
    rng = np.random.default_rng(seed)

    def phi(w):
        Wm = w.reshape(6, 6)
        tot = 0.0
        g = np.zeros((6, 6))
        for key, (its, Vp) in projs.items():
            coef = np.array([Wm[i, a] * Wm[j, b] + (Wm[i, b] * Wm[j, a] if a != b else 0.0) * (1.0) for (i, j, a, b, P) in its])
            # l_i l_j = sum_{a,b} w_ia w_jb x_a x_b ; monomial x_a x_b (a<=b) coefficient = w_ia w_jb + w_ib w_ja (a<b) or w_ia w_ja (a==b)
            # and for i<j the pair appears once (F_ij = F_ji counted once with factor 2 in sum_{i,j}); use factor 2 for i<j
            fac = np.array([2.0 if i != j else 1.0 for (i, j, a, b, P) in its])
            v = Vp @ (coef * fac)
            tot += float(v @ v)
            # gradient
            gv = 2.0 * (Vp.T @ v) * fac
            for t, (i, j, a, b, P) in enumerate(its):
                g[i, a] += gv[t] * Wm[j, b]; g[j, b] += gv[t] * Wm[i, a]
                if a != b:
                    g[i, b] += gv[t] * Wm[j, a]; g[j, a] += gv[t] * Wm[i, b]
        return tot, g.ravel()

    def fg(z):
        n = np.linalg.norm(z); w = z / n
        f, g = phi(w)
        g = (g - (g @ w) * w) / n
        return f, g
    best = []
    for it in range(n_restarts):
        z0 = rng.standard_normal(36); z0 /= np.linalg.norm(z0)
        res = minimize(fg, z0, jac=True, method='L-BFGS-B', options={'maxiter': 500, 'gtol': 1e-14, 'ftol': 1e-18})
        w = res.x / np.linalg.norm(res.x)
        f, _ = phi(w)
        # is D = sum l_i F_i divisible by F (spurious: l_i = c x_i gives D = 6cF)?  check projection of w onto the Euler direction
        eul = np.eye(6).ravel(); eul /= np.linalg.norm(eul)
        best.append((f, float(abs(w @ eul)), w))
    best.sort(key=lambda t: t[0])
    fs = np.array([b[0] for b in best])
    return {'p': p, 'n_items': len(items), 'sectors': len(bykey), 'total_class_rank': tot_rank, 'dimR18_touched': dimR,
            'rank_info': rank_info, 'numeric': {'restarts': n_restarts, 'min_phi': float(fs.min()), 'median_phi': float(np.median(fs)),
            'n_below_1e-8': int(np.sum(fs < 1e-8)), 'best_euler_overlap': best[0][1], 'best_w': [round(float(x), 6) for x in best[0][2]]},
            'seconds': round(time.time() - t0, 1)}


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'd5'
    if which == 'd5':
        for p in (W.PRIME, core.PRIMES[1]):
            t0 = time.time()
            r = d5_test(p)
            r['seconds'] = round(time.time() - t0, 1)
            print(json.dumps(r, indent=None), flush=True)
            W.log(step='hessian-relaxation-d5', p=p, independent=r['independent'], total_rank=r['total_rank'], n_classes=r['n_classes'],
                  zero_products=r['zero_products'], sectors=r['sectors'], seconds=r['seconds'],
                  note=('d=5 EXCLUDED rigorously: the 21 classes [A F_ij] are independent in R_16 mod p, so w_i w_j = 0 for all i,j and w = 0'
                        if r['independent'] else 'd=5 Hessian relaxation has a nontrivial solution space; full E2 needed'))
    elif which == 'd6':
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 200
        r = d6_jtype(W.PRIME, n_restarts=n)
        print(json.dumps({k: v for k, v in r.items() if k != 'rank_info'}, indent=None), flush=True)
        print('rank_info:', [(x['key'], x['dimR'], x['n_classes'], x['rank_classes']) for x in r['rank_info']], flush=True)
        W.log(step='hessian-relaxation-d6-jtype', **{k: v for k, v in r.items() if k != 'rank_info'},
              note='J-type d=6 (D = sum l_i F_i, 36 params, mixed sectors): rank of the class map and numeric minimization of |Phi(w)|^2 over the sphere')
