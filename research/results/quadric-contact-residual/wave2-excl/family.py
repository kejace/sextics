"""E1-admissible families per character sector with explicit gauge, and their Hessian-relaxation closure."""
import sys, time, json
from w2lib import *

def gauge_and_jacobian_family(d, c):
    """Basis of {D in S_d[c] : A D in J} = J_d[c] (+ extras), ordered: gauge elements s*F (s in S_{d-6}[c]) first,
    then independent m*F_i, then non-Jacobian extras.  Returns (Ds, Ns, n_gauge, n_jac, n_extra)."""
    Sd = sector(d, c)
    INV6 = pow(6, p - 2, p)
    Ds, Ns, vecs = [], [], []
    def try_add(D, N):
        v = Sd.vec(D)
        test = nmod_mat(len(vecs) + 1, Sd.n, [x for u in vecs + [v] for x in u], p)
        if test.rank() == len(vecs) + 1:
            vecs.append(v); Ds.append(D); Ns.append(N); return True
        return False
    n_gauge = 0
    if d >= 6:
        for s in monomials(d - 6):
            if chi(s) != c: continue
            sP = {s: 1}
            N = [scale(mul(sP, {tuple(1 if j == i else 0 for j in range(NV)): 1}, p), INV6, p) for i in range(NV)]
            N = [mul(Ap, x, p) for x in N]
            if try_add(mul(sP, Fp, p), N): n_gauge += 1
    n_jac = 0
    if d >= 5:
        for mm in monomials(d - 5):
            for i in range(NV):
                if chi_add(chi(mm), chi_neg(CHI_Z[i])) != c: continue
                N = [{} for _ in range(NV)]; N[i] = mul(Ap, {mm: 1}, p)
                if try_add(shift(Fi[i], mm, 1, p), N): n_jac += 1
    assert len(vecs) == Sd.rank, (len(vecs), Sd.rank)
    # extras
    prof, extras = e1_profile_sector(d, c)
    n_extra = 0
    S6 = sector(d + 6, chi_add(c, CHI_A))
    for P in extras:
        N = S6.lift_field(mul(Ap, P, p))
        assert N is not None
        if try_add(P, N): n_extra += 1
    return Ds, Ns, n_gauge, n_jac, n_extra

def e1_profile_sector(d, c):
    mons = [e for e in monomials(d) if chi(e) == c]
    S6 = sector(d + 6, chi_add(c, CHI_A))
    NF = S6.nf_matrix([S6.vec(shift(Ap, e, 1, p)) for e in mons])
    X, nul = NF.transpose().nullspace()
    Sd = sector(d, c)
    extra = nul - Sd.rank
    chosen = []
    if extra > 0:
        kers = [{e: int(X[i, j]) for i, e in enumerate(mons) if int(X[i, j])} for j in range(nul)]
        NFk = Sd.nf_matrix([Sd.vec(P) for P in kers])
        cur = []
        for i, P in enumerate(kers):
            test = nmod_mat(len(cur) + 1, NFk.ncols(), [int(NFk[t, j]) for t in cur + [i] for j in range(NFk.ncols())], p)
            if test.rank() == len(cur) + 1:
                cur.append(i); chosen.append(P)
            if len(cur) == extra: break
    return {'dimS': len(mons), 'dimJ': Sd.rank, 'kernel': nul, 'extra': extra}, chosen

def sector_closure(d, c, prefix, run_singular=True, timeout=3600, verbose=True):
    """Gauge-fixed Hessian-relaxation closure of the full E1-admissible family in sector c at degree d."""
    t0 = time.time()
    Ds, Ns, ng, nj, ne = gauge_and_jacobian_family(d, c)
    n = len(Ds)
    info = {'d': d, 'chi': str(c), 'n_params': n, 'n_gauge': ng, 'n_jac': nj, 'n_extra': ne}
    if n - ng <= 0:
        info.update({'result': 'no parameters beyond gauge', 'seconds': round(time.time() - t0, 1)})
        return info
    quads, qi = quadrics_from_family(Ds, Ns, d, verbose=verbose)
    # restrict to non-gauge unknowns (gauge coefficients set to 0 -- cross terms with gauge vanish mod J anyway)
    keep = list(range(ng, n))
    quads2 = []
    for cf in quads:
        cf2 = {(k - ng, l - ng): v for (k, l), v in cf.items() if k >= ng and l >= ng}
        if cf2: quads2.append(cf2)
    m = n - ng
    rank, forced, pairs = quad_rank_and_squares(quads2, m)
    info.update({'n_unknowns': m, 'n_quadrics': len(quads2), 'n_target_coords': qi['n_target_coords'], 'quadric_rank': rank,
                 'n_pair_monomials': len(pairs), 'forced_zero_by_squares': forced})
    neq, nv = write_systems(quads2, m, prefix)
    if run_singular and neq > 0:
        v = singular_verdict(prefix + '.sing', timeout=timeout)
        info['singular'] = {k: v[k] for k in ('sdim', 'svdim', 'rc', 'seconds')}
        info['result'] = ('EXCLUDED (only c=0)' if v['sdim'] == 0 else ('positive-dimensional' if (v['sdim'] or 0) > 0 else 'singular-failed'))
    elif neq == 0:
        info['result'] = 'NO EQUATIONS'
    info['seconds'] = round(time.time() - t0, 1)
    with open(prefix + '.family.json', 'w') as f:
        json.dump({'info': info, 'D': [core.to_str(D) for D in Ds]}, f)
    return info
