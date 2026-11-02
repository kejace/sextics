"""Wave-2 task (4): structured irreducible-type denominator candidates (5-variable, Y2-free) at degrees 6..12.

Per candidate: parity check, E1 (A D in J_f, mod p), served-target analysis (dimension of served degree-6 targets
modulo the trivially served Jacobian part), exact confirmation when A is served.  Appends to wave2-search.log.jsonl
and writes <job>.struct.jsonl.  Skips candidates already present in the output file (resumable).
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '2')
import sys, json, time, argparse, itertools, random
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew')
import numpy as np
from fractions import Fraction as Fr
import w2lib as W
from w2lib import core, mul, add, sub, scale, pderiv

P = core.parse


def candidates(prob):
    f = prob.f
    X = [P(v) for v in core.VARS[:5]]
    inv = {
        'G0': P('X0**6+X1**6+X2**6'), 'c': P('X0**3*X1+X1**3*X2'), 'e': P('X0*X1'), 'p': P('Y0*Y1'),
        't': P('Y0**2-Y1**2'), 'm': P('Y0**2+Y1**2'), 's': P('Y0**3-Y1**3'), 'S': P('Y0**3+Y1**3'),
        'Y6': P('Y0**6+Y1**6'), 'f': f, 'cp': P('(X0**3*X1+X1**3*X2)*Y0*Y1'), 'skew': P('X0*X1*Y0*Y1*(Y0**2-Y1**2)/7'),
        'c1': P('X0**3*X1*Y0*Y1'), 'c2': P('X1**3*X2*Y0*Y1'), 'X0': X[0], 'X1': X[1], 'X2': X[2], 'Y0': X[3], 'Y1': X[4],
        'q02': P('X0**2+X2**2'), 'q0y': P('X0**2+Y0**2'), 'x3': P('X0**3+X2**3'), 'y03': P('Y0**3'), 'x03': P('X0**3'),
        'u': P('X0**2'), 'v': P('Y0**2'), 'w': P('X2**2'), 'n': P('X1**2+Y1**2'), 'k': P('X1*Y1'),
    }
    G0, c, e, p, t, m, s, S, Y6, cp, skew = (inv[k] for k in ('G0', 'c', 'e', 'p', 't', 'm', 's', 'S', 'Y6', 'cp', 'skew'))
    out = []

    def addc(name, poly, fam):
        if poly:
            out.append((name, poly, fam))
    # companion pieces of the accepted one-reflection certificate (fqlx) and skew-twisted analogues
    fc = add(add(G0, mul(c, p)), scale(mul(p, mul(p, p)), 2))
    h = add(add(scale(mul(c, c), -1), scale(mul(G0, p), 18)), scale(mul(c, mul(p, p)), 6))
    delta = add(scale(mul(G0, G0), 324), scale(mul(c, mul(c, c)), 24))
    addc('fc', fc, 'companion:fc'); addc('h', h, 'companion:h'); addc('delta', delta, 'companion:delta')
    addc('h+skew*e', add(h, mul(skew, e)), 'companion:h-skewed'); addc('h-skew*e', sub(h, mul(skew, e)), 'companion:h-skewed')
    addc('fc+skew', add(fc, skew), 'companion:fc-skewed')
    addc('f*h/?: h*p', mul(h, p), 'companion:h-product'); addc('h*t', mul(h, t), 'companion:h-product')
    addc('h*m', mul(h, m), 'companion:h-product'); addc('h*e', mul(h, e), 'companion:h-product')
    addc('h*c', mul(h, c), 'companion:h-product'); addc('h*G0', mul(h, G0), 'companion:h-product')
    addc('h*h', mul(h, h), 'companion:h-product')
    addc('delta+h*h*?: delta-324*fc*fc', sub(delta, scale(mul(fc, fc), 324)), 'companion:delta-variant')
    # binary-cubic covariants of P(U) = G0 U^3 + cp U^2 + skew U + Y6 homogenised (degree 12)
    hess = add(add(add(mul(cp, cp), mul(cp, skew)), mul(skew, skew)), add(add(scale(mul(G0, skew), -3), scale(mul(G0, Y6), -9)), scale(mul(cp, Y6), -3)))
    addc('hessP', hess, 'binary-cubic-covariant')
    addc('hessP:b2-3ac', sub(mul(cp, cp), scale(mul(G0, skew), 3)), 'binary-cubic-covariant')
    addc('hessP:bc-9ad', sub(mul(cp, skew), scale(mul(G0, Y6), 9)), 'binary-cubic-covariant')
    addc('hessP:c2-3bd', sub(mul(skew, skew), scale(mul(cp, Y6), 3)), 'binary-cubic-covariant')
    addc('hessP:b2-3ac+c2-3bd', add(sub(mul(cp, cp), scale(mul(G0, skew), 3)), sub(mul(skew, skew), scale(mul(cp, Y6), 3))), 'binary-cubic-covariant')
    # 3x3 principal Hessian minors of f (degree 12) and a few mixed ones
    H = [[pderiv(prob.fi[i], j) for j in range(5)] for i in range(5)]

    def det3(I, J):
        a = [[H[I[r]][J[cc]] for cc in range(3)] for r in range(3)]
        d_ = {}
        for perm, sg in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
            term = mul(mul(a[0][perm[0]], a[1][perm[1]]), a[2][perm[2]])
            d_ = add(d_, scale(term, sg))
        return d_
    for I in itertools.combinations(range(5), 3):
        addc(f'Hess3[{I}]', det3(I, I), 'struct:hessian-3x3')
    for I, J in [((0, 1, 3), (0, 1, 4)), ((0, 2, 3), (0, 2, 4)), ((0, 1, 2), (0, 3, 4)), ((1, 2, 3), (1, 2, 4))]:
        addc(f'Hess3[{I},{J}]', det3(I, J), 'struct:hessian-3x3-mixed')
    # 2x2 Hessian minors times invariants (degrees 10, 12)
    for i, j in [(0, 3), (0, 1), (1, 4), (3, 4), (0, 2), (2, 3)]:
        mnr = sub(mul(H[i][i], H[j][j]), mul(H[i][j], H[i][j]))
        for nm in ('p', 't', 'm', 'e', 'q02', 'n'):
            addc(f'Hess2[{i},{j}]*{nm}', mul(mnr, inv[nm]), 'struct:hessian-2x2-product')
        addc(f'Hess2[{i},{j}]*{"c"}', mul(mnr, c), 'struct:hessian-2x2-product')
    # block radial derivatives times invariants, and partial products at higher degree
    E_blocks = {'rX': add(add(mul(X[0], prob.fi[0]), mul(X[1], prob.fi[1])), mul(X[2], prob.fi[2])),
                'rY': add(mul(X[3], prob.fi[3]), mul(X[4], prob.fi[4])),
                'r0': add(mul(X[0], prob.fi[0]), mul(X[3], prob.fi[3])), 'r14': add(mul(X[1], prob.fi[1]), mul(X[4], prob.fi[4]))}
    for bn, bp in E_blocks.items():
        for nm in ('p', 't', 'm', 'e', 'q02', 'u', 'v', 'w'):
            addc(f'{bn}*{nm}', mul(bp, inv[nm]), 'struct:radial-product')
        for nm in ('G0', 'Y6', 'cp', 'skew', 'f'):
            addc(f'{bn}*{nm}', mul(bp, inv[nm]), 'struct:radial-product')
    for i in range(5):
        for j in range(i, 5):
            pp = mul(prob.fi[i], prob.fi[j])
            for nm in ('p', 't', 'm', 'e', 'u', 'v'):
                addc(f'f_{i}*f_{j}*{nm}', mul(pp, inv[nm]), 'struct:partial-product-quadric')
    # products of invariants at degrees 9..12 (k = 2, 3, 4 factors), rho-even
    names = ['G0', 'c', 'e', 'p', 't', 'm', 's', 'S', 'Y6', 'f', 'cp', 'skew', 'q02', 'q0y', 'x3', 'u', 'v', 'w', 'n', 'k']
    seen = set()
    for k in (2, 3, 4):
        for combo in itertools.combinations_with_replacement(names, k):
            Pp = {(0,) * 6: Fr(1)}
            for nm in combo:
                Pp = mul(Pp, inv[nm])
            d = core.degree(Pp)
            if d < 9 or d > 12:
                continue
            key = tuple(sorted(Pp.items()))
            if key in seen:
                continue
            seen.add(key)
            addc('*'.join(combo), Pp, f'invariant-product:k{k}')
    # cone-type sums (norm-like): f*q + lam*g^2 with g quartic/quintic from the invariant list
    for lam in (Fr(1), Fr(-1), Fr(1, 3), Fr(-1, 3), Fr(3), Fr(-3)):
        for gn in ('u*u', 'v*v', 'p*p', 't*t', 'm*m', 'e*e', 'c', 'cp', 'skew'):
            g = {(0,) * 6: Fr(1)}
            for nm in gn.split('*'):
                g = mul(g, inv[nm])
            dg = core.degree(g)
            if dg == 4:
                addc(f'f*u+({lam})*({gn})^2', add(mul(f, inv['u']), scale(mul(g, g), lam)), 'cone-norm:fu+lam*g^2')
                addc(f'f*v+({lam})*({gn})^2', add(mul(f, inv['v']), scale(mul(g, g), lam)), 'cone-norm:fv+lam*g^2')
            elif dg == 6:
                addc(f'f*f+({lam})*({gn})^2', add(mul(f, f), scale(mul(g, g), lam)), 'cone-norm:f^2+lam*g^2')
    # Fermat-type cones over curves (3-variable and 4-variable), degrees 6..12, rho-even
    for d in (6, 8, 9, 10, 12):
        for keep in itertools.combinations(range(5), 3):
            terms = {}
            for i in keep:
                ex = [0] * 6; ex[i] = d
                if W.rpar(tuple(ex)) == 0:
                    terms[tuple(ex)] = Fr(1)
            if len(terms) >= 2:
                addc(f'cone3-{"".join(core.VARS[i] for i in keep)}-fermat{d}', terms, 'cone:3var-fermat')
                # twisted: add a mixed monomial
                i0, i1 = keep[0], keep[1]
                ex = [0] * 6; ex[i0] = d - 2; ex[i1] = 2
                if W.rpar(tuple(ex)) == 0:
                    tw = dict(terms); tw[tuple(ex)] = Fr(1)
                    addc(f'cone3-{"".join(core.VARS[i] for i in keep)}-twist{d}', tw, 'cone:3var-twisted')
    # random generic calibration points at d = 9, 10 (expected: served_mod_J = 0)
    rng = random.Random(7)
    for d in (9, 10):
        Pp = {mm: Fr(rng.randint(-3, 3)) for mm in prob.Dmons(d, 0, 0)}
        addc(f'random-dense-d{d}', {e_: v for e_, v in Pp.items() if v}, 'random-generic')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--job', default='w2-struct')
    ap.add_argument('--dmin', type=int, default=6)
    ap.add_argument('--dmax', type=int, default=12)
    ap.add_argument('--families', default='', help='comma list of family prefixes to include')
    ap.add_argument('--max-per-degree', type=int, default=0)
    ap.add_argument('--inst', default='one-reflection-skew')
    a = ap.parse_args()
    prob = W.Problem.from_instance(a.inst)
    C = candidates(prob)
    fams = [x for x in a.families.split(',') if x]
    C = [(n, Pp, fam) for n, Pp, fam in C if a.dmin <= core.degree(Pp) <= a.dmax and (not fams or any(fam.startswith(x) for x in fams))]
    C.sort(key=lambda x: core.degree(x[1]))
    if a.max_per_degree:
        byd = {}
        C2 = []
        for n, Pp, fam in C:
            d = core.degree(Pp)
            if byd.get(d, 0) < a.max_per_degree:
                C2.append((n, Pp, fam)); byd[d] = byd.get(d, 0) + 1
        C = C2
    out_path = f'{W.WD}/{a.job}.struct.jsonl'
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            try:
                done.add(json.loads(line)['name'])
            except Exception:
                pass
    print(f'{len(C)} candidates ({len(done)} done)', flush=True)
    W.log_event({'event': 'struct-start', 'job': a.job, 'n_candidates': len(C), 'dmin': a.dmin, 'dmax': a.dmax})
    for name, Pp, fam in C:
        if name in done:
            continue
        t0 = time.time()
        d = core.degree(Pp)
        rec = {'name': name, 'family': fam, 'd': d, 'job': a.job, 'event': 'struct-candidate', 'inst': a.inst}
        pD, sD = W.parity_of(Pp, W.rpar), (W.parity_of(Pp, W.spar) if prob.use_sigma else None)
        rec['rho'], rec['sigma'] = pD, sD
        try:
            if pD is None:
                rec['note'] = 'mixed rho-parity; skipped'
            else:
                rec['e1'] = W.e1_exact_modp(prob, Pp)
                if rec['e1']:
                    sa = W.served_analysis(prob, Pp, pD, sD)
                    rec.update({k: sa.get(k) for k in ('served_dim', 'served_mod_J', 'serves_A_mod_J', 'rank_N', 'nullity_N', 'rows', 'cols')})
                    if sa.get('serves_A_mod_J'):
                        st, inf = W.exact_confirm(prob, Pp, pD, sD, prob.A)
                        rec['confirm'] = st
                        if st == 'solved':
                            ok = W.save_witness_candidate(prob, inf['D6'], inf['N6'], {}, extra={'job': a.job, 'D5': core.to_str(Pp), 'name': name})
                            rec['verified6'] = ok
                            W.log_event({'event': 'WITNESS', 'job': a.job, 'd': d, 'D5': core.to_str(Pp), 'verified6': ok, 'name': name})
                            print('!!!!!!!! WITNESS FOUND !!!!!!!!', name, flush=True)
                else:
                    rec['note'] = 'e1 fails (A D not in J_f mod p)'
        except Exception as ex:
            rec['error'] = f'{ex.__class__.__name__}: {ex}'
        rec['secs'] = round(time.time() - t0, 1)
        rec['D'] = core.to_str(Pp)[:300]
        W.log_event(rec)
        with open(out_path, 'a') as fh:
            fh.write(json.dumps(rec, default=str) + '\n')
        print(f"{name[:50]:50s} d={d} ({pD},{sD}) e1={rec.get('e1')} served={rec.get('served_dim')}/{rec.get('served_mod_J')} "
              f"A={rec.get('serves_A_mod_J')} {rec['secs']}s {rec.get('note', '')}{rec.get('error', '')}", flush=True)
    W.log_event({'event': 'struct-end', 'job': a.job})
    print('DONE', flush=True)


if __name__ == '__main__':
    main()
