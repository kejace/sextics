#!/usr/bin/env python3
"""Mixed-parity Hessian relaxation at degree 6 (aggregator cross-check of the lens' "pure-parity only" gap).

For a lane whose E1 kernel mod J vanishes in degree 6 for BOTH parities (so every E1-admissible degree-6 D is
D = w(F) with w a linear vector field, 36 coefficients c_ik: w_i = sum_k c_ik z_k), the E2-mod-J necessary condition
    A * w^T Hess(F) w in J_12        <=>  (by injectivity of [A] on R_6^even and R_6^odd)   w^T Hess(F) w in J_6
is a system of quadrics in the 36 unknowns c_ik (both parity parts of w^T H w must lie in J_6 separately, since J is
tau-graded).  The Euler direction w = E (c_ik = delta_ik) is the illegal gauge (D = 6F); it is removed by the slice
c_55 = 0 (coefficient of Y2 in w_{Y2}).  "Only the trivial solution" on the slice, certified mod p by an m-primary
leading ideal (msolve -g 1: a pure power of every variable), transfers to C when dim J_6 mod p is the c.i. value
(checked here).  Output: an .ms file (powers written as ^2 -- never c*c) and a JSON summary.

usage: python3 hessrelax6.py <lane> [p]     -> writes hessrelax6-<lane>.ms / .json in the cwd
       python3 hessrelax6.py --check <msolve-output> <n_vars>   -> prints whether every variable has a pure power among the leading monomials
"""
import sys, json, itertools, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')

def build(lane, p=1000000007):
    from hodge import core, jacobian
    from hodge.core import NV
    F, A, q, rec = core.load_instance(lane)
    Fp = core.reduce_mod(F, p)
    # E1 kernel check in degree 6 (both parities): must be 0 for the reduction "D = w(F)" to be complete
    kinfo = {}
    for par in (0, 1):
        out = jacobian.mult_kernel(F, A, 6, p, parity=par)
        kinfo[par] = {'dim': out['dim'], 'ring_dim': out['ring_dim'], 'jacobian_dim': out['jacobian_dim']}
    # reducers in degree 6 (both parities); c.i. values: R_6 = 218 (even) + 208 (odd), J_6 = 18 + 18
    red = {par: jacobian._Reducer(F, 6, p, par) for par in (0, 1)}
    hilbert_ok = (red[0].ring_dim == 218 and red[1].ring_dim == 208)
    # second partials
    Fi = [core.pderiv(Fp, i, p) for i in range(NV)]
    Fij = [[core.pderiv(Fi[i], j, p) for j in range(NV)] for i in range(NV)]
    # variables c_{ik}: index v = 6*i + k ; w_i = sum_k c_ik z_k
    # w^T H w = sum_{i,j} w_i w_j F_ij = sum_{i,k,j,l} c_ik c_jl z_k z_l F_ij
    # quadric coefficients: for each pair-monomial c_ik c_jl (unordered), poly = z_k z_l F_ij (+ z_l z_k F_ji if different index pair)
    nvar = NV * NV
    pairs = {}  # (v1<=v2) -> poly (mod p)
    for i in range(NV):
        for k in range(NV):
            for j in range(NV):
                for l in range(NV):
                    v1, v2 = NV * i + k, NV * j + l
                    key = (min(v1, v2), max(v1, v2))
                    e = [0] * NV; e[k] += 1; e[l] += 1
                    poly = core.shift(Fij[i][j], tuple(e), 1, p)
                    pairs[key] = core.add(pairs.get(key, {}), poly, p)
    # normal forms and quadric extraction: for each parity, coordinates = standard monomials
    quadrics = []  # list of dict {pair-key: coeff}
    nf_cache = {}
    for par in (0, 1):
        std = red[par].standard
        sidx = {m: t for t, m in enumerate(std)}
        rows = [dict() for _ in std]
        for key, poly in pairs.items():
            part = {e: c for e, c in poly.items() if core.yparity(e) == par}
            if not part:
                continue
            nf = red[par].normal_form(part)
            for e, c in nf.items():
                rows[sidx[e]][key] = c
        quadrics.extend([r for r in rows if r])
    # gauge slice: c_55 = 0  (variable index 35); drop that variable from all quadrics
    gauge_var = NV * 5 + 5
    names = [f'c{i}_{k}' for i in range(NV) for k in range(NV)]
    keep = [v for v in range(nvar) if v != gauge_var]
    def term(key, c):
        v1, v2 = key
        if v1 == gauge_var or v2 == gauge_var:
            return None
        if v1 == v2:
            return f'{c}*{names[v1]}^2'
        return f'{c}*{names[v1]}*{names[v2]}'
    polys = []
    for qd in quadrics:
        terms = [t for t in (term(k, c) for k, c in qd.items()) if t]
        if terms:
            polys.append('+'.join(terms))
    polys = sorted(set(polys))
    ms = ','.join(names[v] for v in keep) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n'
    open(f'hessrelax6-{lane}.ms', 'w').write(ms)
    summary = {'lane': lane, 'p': p, 'e1_kernel_d6': kinfo, 'hilbert_ok': hilbert_ok, 'ring_dims': {par: red[par].ring_dim for par in (0, 1)},
               'n_vars': len(keep), 'n_quadrics': len(polys), 'gauge': 'c5_5 = 0 (Euler direction removed)', 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}
    json.dump(summary, open(f'hessrelax6-{lane}.json', 'w'), indent=1)
    print(json.dumps(summary))

def check(path, names):
    txt = open(path).read()
    # msolve -g 1 output: list of leading monomials like [c0_0^2, c0_1*c0_2, ...]
    import re
    body = txt[txt.find('['):]
    monos = [m.strip() for m in body.strip().strip('[]:').split(',') if m.strip()]
    pure = {}
    for m in monos:
        m2 = m.replace(' ', '')
        mm = re.fullmatch(r'(?:1\*)?([A-Za-z0-9_]+)(?:\^(\d+))?', m2)
        if mm:
            pure[mm.group(1)] = int(mm.group(2) or 1)
    missing = [n for n in names if n not in pure]
    print(json.dumps({'n_leading': len(monos), 'pure_powers': len(pure), 'missing_pure_power_for': missing, 'm_primary': not missing,
                      'max_pure_power': max(pure.values()) if pure else None}))
    return not missing

if __name__ == '__main__':
    if sys.argv[1] == '--check':
        names = [f'c{i}_{k}' for i in range(6) for k in range(6)]
        names.remove('c5_5')
        check(sys.argv[2], names)
    else:
        build(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000000007)
