"""w2s_structured.py -- structured irreducible-type candidates at degrees 10 and 12 (Y-even):
numerical E1 screen (projection of A*D onto the complement of J_{d+6}^-), exact mod-p E1 for the survivors,
then the full normalized-E2 residual (NormE2 with the single-D family) and exact fixed-D solve if it is ~0."""
import sys, os, time, json, itertools
from functools import reduce
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from w2s_common import *
from w2s_nvp import NormE2


def blocks(F):
    out = {0: {}, 2: {}, 4: {}, 6: {}}
    for e, c in F.items():
        out[e[3] + e[4] + e[5]][e] = c
    return out


def hess_det(P, idx):
    """determinant of the Hessian of P restricted to the index list (exact, Fractions)."""
    H = [[pderiv(pderiv(P, i), j) for j in idx] for i in idx]
    n = len(idx)
    det = {}
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for a in range(n):
            for b in range(a + 1, n):
                if perm[a] > perm[b]: sgn = -sgn
        term = {(0,) * NV: Fraction(sgn)}
        for a in range(n):
            term = mul(term, H[a][perm[a]])
            if not term: break
        det = add(det, term)
    return det


def candidates(F, A):
    B = blocks(F)
    a, b, c, d_ = B[0], B[2], B[4], B[6]
    Fi = [pderiv(F, i) for i in range(NV)]
    Z = [core.parse(v) for v in core.VARS]
    C = []
    # degree 12: Hessian covariants of the binary cubic a + b t + c t^2 + d t^3 of the cubic cover X/tau -> P2 x P2
    C.append(('d12 cubic-cover Hessian covariant b^2-3ac', sub(mul(b, b), scale(mul(a, c), 3))))
    C.append(('d12 cubic-cover Hessian covariant bc-9ad', sub(mul(b, c), scale(mul(a, d_), 9))))
    C.append(('d12 cubic-cover Hessian covariant c^2-3bd', sub(mul(c, c), scale(mul(b, d_), 3))))
    C.append(('d12 cubic-cover H = (b^2-3ac)+(bc-9ad)+(c^2-3bd)', add(add(sub(mul(b, b), scale(mul(a, c), 3)), sub(mul(b, c), scale(mul(a, d_), 9))), sub(mul(c, c), scale(mul(b, d_), 3)))))
    C.append(('d12 b^2-4ac', sub(mul(b, b), scale(mul(a, c), 4))))
    C.append(('d12 c^2-4bd', sub(mul(c, c), scale(mul(b, d_), 4))))
    C.append(('d12 det Hess_X(a)', hess_det(a, [0, 1, 2])))
    C.append(('d12 det Hess_Y(d)', hess_det(d_, [3, 4, 5])))
    for trip in itertools.combinations(range(NV), 3):
        C.append((f'd12 3x3 principal Hessian minor {trip}', hess_det(F, list(trip))))
    # degree 12: sum of squared Euler cross terms, sum over pairs (F_i Z_j - F_j Z_i)^2
    S = {}
    for i in range(NV):
        for j in range(i + 1, NV):
            t = sub(mul(Fi[i], Z[j]), mul(Fi[j], Z[i])); S = add(S, mul(t, t))
    C.append(('d12 sum_{i<j} (F_i Z_j - F_j Z_i)^2', S))
    # degree 10: gradient norms (Jacobian-type, irreducible)
    C.append(('d10 sum F_i^2', reduce(add, (mul(Fi[i], Fi[i]) for i in range(NV)))))
    C.append(('d10 sum_X F_i^2 - sum_Y F_i^2', sub(reduce(add, (mul(Fi[i], Fi[i]) for i in range(3))), reduce(add, (mul(Fi[i], Fi[i]) for i in range(3, 6))))))
    C.append(('d10 F_X0^2 + F_Y0^2', add(mul(Fi[0], Fi[0]), mul(Fi[3], Fi[3]))))
    C.append(('d10 F_X0^2 - F_Y0^2', sub(mul(Fi[0], Fi[0]), mul(Fi[3], Fi[3]))))
    C.append(('d10 sum Z_i^5 F_i', reduce(add, (mul(power(Z[i], 5), Fi[i]) for i in range(NV)))))
    # degree 10: second polars w^T Hess(F) w for structured cubic fields w
    H = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
    def wHw(w):
        out = {}
        for i in range(NV):
            if not w[i]: continue
            for j in range(NV):
                if not w[j]: continue
                out = add(out, mul(mul(w[i], w[j]), H[i][j]))
        return out
    C.append(('d10 wHw, w = grad(sum Z_i^4)', wHw([scale(power(Z[i], 3), 4) for i in range(NV)])))
    C.append(('d10 wHw, w = grad(X0^4 + Y0^4)', wHw([scale(power(Z[0], 3), 4), {}, {}, scale(power(Z[3], 3), 4), {}, {}])))
    C.append(('d10 wHw, w = grad(X0^2 Y0^2)', wHw([mul(scale(Z[0], 2), power(Z[3], 2)), {}, {}, mul(scale(Z[3], 2), power(Z[0], 2)), {}, {}])))
    C.append(('d10 wHw, w = J grad(sum Z_i^4) (symplectic pairs)', wHw([scale(power(Z[3], 3), 4), scale(power(Z[4], 3), 4), scale(power(Z[5], 3), 4), scale(power(Z[0], 3), -4), scale(power(Z[1], 3), -4), scale(power(Z[2], 3), -4)])))
    C.append(('d10 wHw, w_i = Z_i^3', wHw([power(Z[i], 3) for i in range(NV)])))
    C.append(('d10 wHw, w = (Y0^3,0,0,X0^3,0,0)', wHw([power(Z[3], 3), {}, {}, power(Z[0], 3), {}, {}])))
    return C


def main():
    F, A, q = load()
    Fi = [vec(pderiv(F, i), 5) for i in range(NV)]
    Av = vec(A, 6)
    out = []
    for name, D in candidates(F, A):
        t0 = time.time()
        if not D:
            continue
        d = core.degree(D); pD = core.parity_of(D)
        rec = dict(name=name, d=d, pD=pD, terms=len(D))
        if pD != 0 or not core.is_homogeneous(D):
            rec['note'] = 'not Y-even homogeneous'; out.append(rec); log_event(event='structured', **rec); continue
        Dv = vec(D, d)
        U, rJ, rows = J_complement(Fi, d + 6, 1, verbose=False)
        AD = polymul(Av, 6, Dv, d)[rows]
        e1num = float(np.linalg.norm(U.T @ AD) / max(np.linalg.norm(AD), 1e-300))
        rec['e1_num_resid'] = e1num
        divF, _ = jacobian.divisible_by(F, D)
        rec['F_divides_D'] = bool(divF)
        if e1num < 1e-8:
            ok, info = fixed_d.e1_solvable(F, A, q, D)
            rec['e1_exact'] = bool(ok); rec['e1_info'] = {k: int(info[k]) for k in ('rank', 'rank_aug')}
            if ok and d <= 12:
                Dm = Dv[None, :]
                N0, res = e1_lifts(F, A, Fi, d, 0, Dm, verbose=False)
                E2 = NormE2(F, A, d, 0, Dm, N0, P=None, label=name, verbose=False)
                f, rel, x, rt = E2.value(np.array([1.0]))
                rec['e2_num_rel'] = float(rel); rec['e1_lift_resid'] = float(res.max())
                if rel < 1e-8:
                    status, data = exact_confirm(F, A, q, D)
                    rec['exact_status'] = status
                    if status == 'solved':
                        w = save_witness_candidate(F, A, q, D, data['N'], data['Q'])
                        log_event(event='WITNESS', name=name, d=d, verified=w['verified'], F_divides_D=w['F_divides_D'])
        rec['seconds'] = round(time.time() - t0, 1)
        out.append(rec)
        log_event(event='structured', **rec)
        print(json.dumps(rec), flush=True)
        with open(os.path.join(LANE, 'w2s_structured.json'), 'w') as f:
            json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
