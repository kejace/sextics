"""E1-admissible denominator spaces per degree and Y-parity (exact mod p): mult_kernel + J(F)_d dims.

  python3 wave2/e1spaces.py --dmin 4 --dmax 10 [--parities 0,1] [--prime-index 0] [--tag e1spaces]
Writes out/<tag>.jsonl (one line per (d, parity)) and out/e1kernel_d{d}_p{par}.json with the kernel basis mod p
(representatives in standard monomials) and the E1 lifts of the kernel representatives.
"""
import sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('--dmin', type=int, default=4)
ap.add_argument('--dmax', type=int, default=10)
ap.add_argument('--parities', default='0,1')
ap.add_argument('--prime-index', type=int, default=0)
ap.add_argument('--tag', default='e1spaces')
ap.add_argument('--no-lifts', action='store_true')
args = ap.parse_args()

F, A, q = load()
p = core.PRIMES[args.prime_index]
Ap = core.reduce_mod(A, p); Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
pA = core.parity_of(A)


def kernel_with_lifts(d, par):
    """{D in S_d^par : A D in J} via nullspace of [A m | J gens]; returns dims, kernel reps mod J, lifts."""
    t0 = time.time()
    pT = (pA + par) % 2
    Dm = monomials(d, par)
    gens = jacobian.jacobian_generators(F, d + 6, p, pT)  # (i, m, poly)
    cols = [shift(Ap, m, 1, p) for m in Dm] + [g[2] for g in gens]
    rows = {}
    for c in cols:
        for e in c:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, c in enumerate(cols):
        for e, v in c.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    K = len(Dm)
    if nullity == 0:
        return dict(d=d, parity=par, rows=R, cols=C, nD=K, nJgens=len(gens), nullity=0, e1_dim=0, J_dim=0, kernel_mod_J=0,
                    t_build=round(t1 - t0, 1), t_null=round(t2 - t1, 1)), [], []
    # project to D coords, rref
    proj = nmod_mat(nullity, K, [int(X[j, k]) for k in range(nullity) for j in range(K)], p)
    Pr, rk = proj.rref()
    # also keep full nullspace vectors to extract lifts: we need, for each kernel rep D, an N with N(F) = A D.
    # Take the rref transform: solve proj^T y = e for each basis vector? Simpler: recompute lift by solving J-gens system later.
    e1_dim = rk
    # J_d dimension (this parity)
    red = jacobian._Reducer(F, d, p, par)
    J_dim = red.jacobian_dim
    # kernel reps modulo J: normal forms of the rref rows, then independent ones
    vecs = []
    for i in range(rk):
        v = {Dm[j]: int(Pr[i, j]) for j in range(K) if int(Pr[i, j])}
        vecs.append(v)
    reduced = [red.normal_form(v) for v in vecs]
    reduced = [v for v in reduced if v]
    basis = []
    if reduced:
        sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
        flat2 = [0] * (len(reduced) * len(sm))
        for r_, v in enumerate(reduced):
            for e, c in v.items():
                flat2[r_ * len(sm) + sidx[e]] = c
        Qm = nmod_mat(len(reduced), len(sm), flat2, p)
        Qr, rk2 = Qm.rref()
        for i in range(rk2):
            basis.append({sm[j]: int(Qr[i, j]) for j in range(len(sm)) if int(Qr[i, j])})
    info = dict(d=d, parity=par, rows=R, cols=C, nD=K, nJgens=len(gens), nullity=int(nullity), e1_dim=int(e1_dim),
                J_dim=int(J_dim), kernel_mod_J=len(basis), ring_dim=red.ring_dim, t_build=round(t1 - t0, 1), t_null=round(t2 - t1, 1))
    lifts = []
    if basis and not args.no_lifts:
        # lifts: solve  sum_g y_g gen_g = A D  for each kernel rep (one augmented rref)
        t3 = time.time()
        rhs = [mul(Ap, b, p) for b in basis]
        rows2 = {}
        for c in [g[2] for g in gens] + rhs:
            for e in c:
                rows2.setdefault(e, len(rows2))
        R2, C2 = len(rows2), len(gens)
        nb = len(rhs)
        flat3 = [0] * (R2 * (C2 + nb))
        for j, g in enumerate(gens):
            for e, v in g[2].items():
                flat3[rows2[e] * (C2 + nb) + j] = v
        for t_, r_ in enumerate(rhs):
            for e, v in r_.items():
                flat3[rows2[e] * (C2 + nb) + C2 + t_] = v
        Maug = nmod_mat(R2, C2 + nb, flat3, p)
        Mr, rank = Maug.rref()
        pivots = []
        for i in range(rank):
            j = 0
            while j < C2 + nb and int(Mr[i, j]) == 0:
                j += 1
            pivots.append(j)
        assert all(pv < C2 for pv in pivots), 'lift inconsistent (should not happen for kernel reps)'
        for t_ in range(nb):
            N = [{} for _ in range(NV)]
            for i, pv in enumerate(pivots):
                v = int(Mr[i, C2 + t_])
                if v:
                    gi, gm, _ = gens[pv]
                    N[gi][gm] = (N[gi].get(gm, 0) + v) % p
            # verify N(F) = A D mod p
            NF = {}
            for i in range(NV):
                NF = add(NF, mul(N[i], Fi[i], p), p)
            assert sub(NF, rhs[t_], p) == {}, 'lift verification failed'
            lifts.append([core.to_str({e: v for e, v in Ni.items()}) for Ni in N])
        info['t_lifts'] = round(time.time() - t3, 1)
    return info, basis, lifts


if __name__ == '__main__':
    pars = [int(x) for x in args.parities.split(',')]
    for d in range(args.dmin, args.dmax + 1):
        for par in pars:
            info, basis, lifts = kernel_with_lifts(d, par)
            info['p'] = p
            info['batch'] = 'e1spaces'
            print(json.dumps(info), flush=True)
            log_line(f'{OUT}/{args.tag}.jsonl', info)
            if basis:
                with open(f'{OUT}/e1kernel_d{d}_p{par}.json', 'w') as fh:
                    json.dump({'info': info, 'p': p, 'basis_modp': [core.to_str(b) for b in basis], 'lifts_modp': lifts}, fh)
