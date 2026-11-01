"""Koszul-parametrized exclusion at degree d (parity-free, (U,V) frame, exact rationals).
D = sum_i a_i P_i  (a_i of degree d-5),  N = A a + sum_{i<j} k_ij (P_j e_i - P_i e_j)  (k_ij of degree d-4),
E2:  D div N - N(D)  ==  0  mod Fhat   (Fhat monic in U0 up to the factor 1/2).
Outputs a Singular script whose ideal I (in the parameters) cuts out all (a,k) with E2; the legal-D question
is whether V(I) contains a point with D != 0 (d=5) / D not proportional to F (d=6)."""
import sys, json, itertools, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from fractions import Fraction as Fr
d = int(sys.argv[1])
uv = json.load(open('uv_forms.json'))
Fh = core.parse(uv['Fhat']); Ah = core.parse(uv['Ahat'])
NV = 6
P = [core.pderiv(Fh, i) for i in range(NV)]
mons_a = core.monomials(d - 5) if d >= 5 else []
mons_k = core.monomials(d - 4) if d >= 4 else []
# parameters: a[i][m], k[(i,j)][m]
params = []
aidx = {}; kidx = {}
for i in range(NV):
    for m in mons_a:
        aidx[(i, m)] = len(params); params.append(f'a{i}_' + ''.join(map(str, m)))
for i in range(NV):
    for j in range(i + 1, NV):
        for m in mons_k:
            kidx[(i, j, m)] = len(params); params.append(f'k{i}{j}_' + ''.join(map(str, m)))
npar = len(params)
print('degree', d, 'parameters', npar, '(a:', len(aidx), ', k:', len(kidx), ')', flush=True)
# polynomials with parameter-linear coefficients: dict {(zmon): {paramindex: coeff}}  (linear in params)
def lin_scale_mon(poly, m, pidx):
    """poly * monomial m, times parameter pidx -> linear form dict"""
    out = {}
    for e, c in poly.items():
        e2 = tuple(e[t] + m[t] for t in range(NV))
        out.setdefault(e2, {})[pidx] = out.get(e2, {}).get(pidx, 0) + c
    return out
def lin_add(A_, B_):
    out = {e: dict(v) for e, v in A_.items()}
    for e, v in B_.items():
        dst = out.setdefault(e, {})
        for pi, c in v.items():
            dst[pi] = dst.get(pi, 0) + c
            if dst[pi] == 0: del dst[pi]
        if not dst: del out[e]
    return out
def lin_deriv(A_, i):
    out = {}
    for e, v in A_.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1; e2 = tuple(e2)
            out[e2] = {pi: c * e[i] for pi, c in v.items()}
    return out
# D = sum_i a_i P_i
Dl = {}
for i in range(NV):
    for m in mons_a:
        Dl = lin_add(Dl, lin_scale_mon(P[i], m, aidx[(i, m)]))
# N_i = A a_i + sum_{j>i} k_ij P_j - sum_{j<i} k_ji P_i ... careful: term k_ij (P_j e_i - P_i e_j)
Nl = [{} for _ in range(NV)]
for i in range(NV):
    for m in mons_a:
        Nl[i] = lin_add(Nl[i], lin_scale_mon(Ah, m, aidx[(i, m)]))
for i in range(NV):
    for j in range(i + 1, NV):
        for m in mons_k:
            pid = kidx[(i, j, m)]
            Nl[i] = lin_add(Nl[i], lin_scale_mon(P[j], m, pid))
            Nl[j] = lin_add(Nl[j], lin_scale_mon(core.scale(P[i], -1), m, pid))
divN = {}
for i in range(NV):
    divN = lin_add(divN, lin_deriv(Nl[i], i))
# quadratic forms: dict {zmon: {(p1,p2) sorted: coeff}}
def quad_mul(A_, B_):
    out = {}
    for e1, v1 in A_.items():
        for e2, v2 in B_.items():
            e = tuple(e1[t] + e2[t] for t in range(NV))
            dst = out.setdefault(e, {})
            for p1, c1 in v1.items():
                for p2, c2 in v2.items():
                    key = (p1, p2) if p1 <= p2 else (p2, p1)
                    dst[key] = dst.get(key, 0) + c1 * c2
    return out
def quad_add(A_, B_, sign=1):
    out = {e: dict(v) for e, v in A_.items()}
    for e, v in B_.items():
        dst = out.setdefault(e, {})
        for key, c in v.items():
            dst[key] = dst.get(key, 0) + sign * c
            if dst[key] == 0: del dst[key]
        if not dst: del out[e]
    return out
t0 = time.time()
E2 = quad_mul(Dl, divN)
for i in range(NV):
    E2 = quad_add(E2, quad_mul(Nl[i], lin_deriv(Dl, i)), -1)
print('E2 built: z-monomials', len(E2), 'in', round(time.time() - t0, 1), 's', flush=True)
# reduce modulo Fhat: U0^6 -> -2*(Fhat - U0^6/2)  (Fhat has U0^6 coefficient 1/2)
assert Fh[(6, 0, 0, 0, 0, 0)] == Fr(1, 2)
rest = core.scale(core.sub(Fh, {(6, 0, 0, 0, 0, 0): Fr(1, 2)}), -2)   # U0^6 == rest  mod Fhat
cache = {}
def reduce_mon(e):
    """normal form of monomial e modulo Fhat (as dict), with U0-degree < 6"""
    if e[0] < 6: return {e: Fr(1)}
    if e in cache: return cache[e]
    e2 = list(e); e2[0] -= 6; e2 = tuple(e2)
    out = {}
    for r_e, r_c in rest.items():
        m = tuple(r_e[t] + e2[t] for t in range(NV))
        for mm, cc in reduce_mon(m).items():
            out[mm] = out.get(mm, 0) + r_c * cc
    out = {k: v for k, v in out.items() if v}
    cache[e] = out
    return out
E2r = {}
for e, v in E2.items():
    for mm, cc in reduce_mon(e).items():
        dst = E2r.setdefault(mm, {})
        for key, c in v.items():
            dst[key] = dst.get(key, 0) + cc * c
            if dst[key] == 0: del dst[key]
        if not dst: del E2r[mm]
print('reduced mod Fhat: equations', len(E2r), 'in', round(time.time() - t0, 1), 's', flush=True)
# write Singular script
def fmt(c):
    c = Fr(c)
    return f'({c.numerator}/{c.denominator})' if c.denominator != 1 else (f'({c.numerator})' if c.numerator < 0 else str(c.numerator))
lines = []
lines.append(f'ring r = 0, ({",".join(params)}), dp;')
lines.append('option(redSB);')
lines.append(f'ideal I;')
eqs = []
for e, v in E2r.items():
    terms = []
    for (p1, p2), c in v.items():
        terms.append(f'{fmt(c)}*{params[p1]}*{params[p2]}')
    eqs.append('+'.join(terms))
for n_, eq in enumerate(eqs):
    lines.append(f'I[{n_+1}] = {eq};')
with open(f'koszul_d{d}.sing', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
json.dump({'params': params, 'npar': npar, 'neq': len(eqs)}, open(f'koszul_d{d}_meta.json', 'w'))
print('wrote', f'koszul_d{d}.sing', 'equations', len(eqs), 'params', npar, flush=True)
