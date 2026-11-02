"""Job H (excl): structured core-only families in the FULL 6-variable system restricted to the forced character
sector (lossless for G-eigenvectors; mod P2, plus exact ranks over Q when D is rational and the restricted system
is small).  Families (all core-only, pure (Y0,Y1)-parity, degree <= 12):
  * cones with vertex at F_p-points of the sextic surface S_H = {H = 0} c P^3 over plane curves (Fermat, random,
    g-shaped) of degrees 2..6 (vertex v: D(x) = C(l1(x), l2(x), l3(x)) with l_i linear forms vanishing at v);
  * cones with vertex a line of P^3 meeting S_H (D in two linear forms) — binary forms in (l1,l2): products of lines, excluded by the sum principle, so only as controls;
  * H + lambda L^6 (six planes through a point of S_H) and other pencils c + lambda H, c + lambda f, c + lambda g with c a kernel rep;
  * powers/products: c*H, c*f etc. are reducible (sum principle) -- skipped;
  * degree-12: f*g, H^2 (reducible), f^2+g^2, (f+g)^2 - H^2, kernel-rep quadratic combos;
  * random core-only polynomials of degrees 6..9 that are NOT E1-admissible (negative controls: E1 must fail).
Each D is first checked for E1 (mod p); E2 only if E1 passes."""
import sys, time, json, random, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log
from hodge import core, symmetry
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV, power, reduce_mod
from fractions import Fraction
from multiprocessing import Pool

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-struct'
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
P = core.parse
p = W.P2
H = W.H
f = {e: c for e, c in H.items() if e[3] == 0 and e[4] == 0}
g = {e: c for e, c in H.items() if e[0] == 0 and e[1] == 0}
K = sub(sub(H, f), g)
Hp = reduce_mod(H, p)

def task(args):
    tag, family, Dj, modp = args
    D = {tuple(int(x) for x in k.split()): (int(v) if modp else Fraction(v)) for k, v in Dj.items()}
    t0 = time.time()
    try:
        d = core.degree(D)
        if core.parity_of(D) is None:
            return {'exp': 'struct-fixedD', 'D': tag, 'family': family, 'd': d, 'prime': p, 'result': 'skipped: mixed parity'}
        ok, info = W.e1_test(D, p)
        if not ok:
            return {'exp': 'struct-fixedD', 'D': tag, 'family': family, 'd': d, 'prime': p, 'e1': False,
                    'result': 'E1 fails (A D not in J(F)) [MODP]', 'seconds': round(time.time() - t0, 1)}
        allowed, sub_ = W.sector_allowed(D)
        S = W.build_fixed_D(D, p, allowed=allowed)
        r = W.modp_consistency(S)
        out = {'exp': 'struct-fixedD', 'D': tag, 'family': family, 'd': d, 'prime': p, 'e1': True, 'stabilizer': len(sub_),
               'rows': r['rows'], 'cols': r['cols'], 'rank_M_modp': r['rank_M'], 'rank_Mb_modp': r['rank_Mb'],
               'consistent_modp': r['consistent'], 'nullity_modp': r['nullity']}
        if r['consistent']:
            out['result'] = 'CONSISTENT mod p -- LIFT NEEDED'
        elif not modp and r['cols'] <= 2500:
            Sx = W.build_fixed_D(D, None, allowed=allowed)
            rM = W.rank_exact(Sx['colpolys']); rMb = W.rank_exact(Sx['colpolys'], extra=[Sx['rhs']])
            out.update({'rank_M_Q': rM, 'rank_Mb_Q': rMb})
            out['result'] = 'inconsistent over C [EXACT]' if rMb > rM else 'CONSISTENT over Q (exact) -- SOLVE NEEDED'
        else:
            out['result'] = 'inconsistent mod p [MODP]'
        out['seconds'] = round(time.time() - t0, 1)
        return out
    except Exception as ex:
        return {'exp': 'struct-fixedD', 'D': tag, 'family': family, 'prime': p, 'result': f'ERROR {ex!r}', 'seconds': round(time.time() - t0, 1)}

def lin_forms_vanishing_at(v):
    """three independent linear forms (mod p) in the core variables vanishing at the point v of P^3 (core coords)."""
    # v = (v0, v1, v3, v4) in F_p^4; forms l(x) = sum c_i x_i with sum c_i v_i = 0: basis of the orthogonal complement
    idx = list(W.CORE)
    piv = next(i for i in range(4) if v[i] % p)
    forms = []
    for i in range(4):
        if i == piv: continue
        # x_i - (v_i / v_piv) x_piv
        c = (-v[i] * pow(v[piv], -1, p)) % p
        l = {tuple(1 if j == idx[i] else 0 for j in range(NV)): 1}
        if c: l[tuple(1 if j == idx[piv] else 0 for j in range(NV))] = c
        forms.append(l)
    return forms

def points_on_SH(n, rng):
    """random F_p-points of S_H: fix X0=1, random X1, Y0, solve the sextic in Y1 by brute force on a random line."""
    pts = []
    tries = 0
    while len(pts) < n and tries < 20000:
        tries += 1
        a, b = rng.randrange(p), rng.randrange(p)
        # restrict to the line (1, a, b, t): H(1,a,b,t) is a univariate sextic in t; find roots via flint
        from flint import nmod_poly
        coeffs = [0] * 7
        for e, c in Hp.items():
            coeffs[e[4]] = (coeffs[e[4]] + c * pow(a, e[1], p) * pow(b, e[3], p)) % p
        poly = nmod_poly(coeffs, p)
        for fac, mult in poly.factor()[1]:
            if fac.degree() == 1:
                t = (-int(fac.coeffs()[0]) * pow(int(fac.coeffs()[1]), -1, p)) % p
                pts.append((1, a, b, t))
                break
    return pts

def compose(C3, forms):
    """C3: poly in 3 'slot' variables encoded as dict {(i,j,k): coeff}; substitute slots by the linear forms."""
    out = {}
    for (i, j, k), c in C3.items():
        term = {(0,) * NV: c % p}
        term = mul(term, power(forms[0], i, p), p) if i else term
        term = mul(term, power(forms[1], j, p), p) if j else term
        term = mul(term, power(forms[2], k, p), p) if k else term
        out = add(out, term, p)
    return out

if __name__ == '__main__':
    t0 = time.time()
    rng = random.Random(99)
    tasks = []
    def addT(tag, family, D, modp):
        tasks.append((tag, family, {' '.join(map(str, e)): str(c) for e, c in D.items()}, modp))
    # kernel reps for pencils
    reps6 = {par: W.small_kernel_reps(W.H, W.CORE, W.A, 6, par)['reps_modJ'] for par in (0, 1)}
    # --- cones with vertex on S_H (mod p)
    pts = points_on_SH(4, rng)
    log(JOB, {'exp': 'struct-points', 'n_points_on_SH': len(pts), 'points': pts, 'prime': p})
    for vi, v in enumerate(pts):
        forms = lin_forms_vanishing_at(v)
        for deg in (2, 3, 4, 5, 6):
            fermat = {(deg, 0, 0): 1, (0, deg, 0): 1, (0, 0, deg): 1}
            addT(f'Fermat cone deg {deg}, vertex #{vi} on S_H', 'cone-vertex-on-SH', compose(fermat, forms), True)
            rand3 = {(i, j, deg - i - j): rng.randrange(1, p) for i in range(deg + 1) for j in range(deg + 1 - i)}
            addT(f'random cone deg {deg}, vertex #{vi} on S_H', 'cone-vertex-on-SH', compose(rand3, forms), True)
        # six planes through the vertex: product of 6 random linear forms through v (reducible: control for sum principle)
        prodl = {(0,) * NV: 1}
        for _ in range(6):
            l = add(add(scale(forms[0], rng.randrange(1, p), p), scale(forms[1], rng.randrange(1, p), p), p), scale(forms[2], rng.randrange(1, p), p), p)
            prodl = mul(prodl, l, p)
        addT(f'six planes through vertex #{vi} on S_H', 'planes-through-point-on-SH', prodl, True)
    # cone with vertex off S_H (control)
    forms = lin_forms_vanishing_at((1, 0, 0, 0))
    addT('Fermat sextic cone, vertex (1:0:0:0) off S_H', 'cone-vertex-off-SH', compose({(6, 0, 0): 1, (0, 6, 0): 1, (0, 0, 6): 1}, forms), True)
    # --- rational structured forms (exact ranks possible)
    addT('f^2 + g^2', 'd12-structural', add(mul(f, f), mul(g, g)), False)
    addT('f*g', 'd12-structural', mul(f, g), False)
    addT('H^2 - f*g', 'd12-structural', sub(mul(H, H), mul(f, g)), False)
    addT('K^2', 'd12-structural', mul(K, K), False)
    addT('f^2 - g^2', 'd12-structural', sub(mul(f, f), mul(g, g)), False)
    addT('X0^6 + X1^6 + Y0^6 + Y1^6 (Fermat core)', 'd6-structural', P('X0**6+X1**6+Y0**6+Y1**6'), False)
    addT('f + g (split core)', 'd6-structural', add(f, g), False)
    addT('H - 2K', 'd6-structural', sub(H, scale(K, 2)), False)
    addT('(X0^2+Y0^2)^3', 'd6-structural', power(P('X0**2+Y0**2'), 3), False)
    addT('X0^2 Y0^2 (X0^2 + 2 Y0^2)', 'd6-structural', P('X0**2*Y0**2*(X0**2+2*Y0**2)'), False)
    addT('X0 X1 Y0 Y1 (X0^2 + 2 Y0^2)', 'd6-structural', P('X0*X1*Y0*Y1*(X0**2+2*Y0**2)'), False)
    # pencils with the first kernel reps
    for par in (0, 1):
        for i, c in enumerate(reps6[par][:3]):
            for name, Q6 in (('H', H), ('f', f), ('g', g), ('K', K)):
                if core.parity_of(add(c, Q6)) is None: continue
                for lam in (1, -1):
                    addT(f'ker6 rep #{i} par{par} + ({lam})*{name}', 'pencil-rep+structural', add(c, scale(Q6, Fraction(lam))), False)
            # degree-12: c*rep' sums (not products): c_i^2 + c_j^2
            for jx, c2 in enumerate(reps6[par][:3]):
                if jx <= i: continue
                addT(f'(ker6 rep #{i})^2 + (ker6 rep #{jx})^2 par{par}', 'd12-rep-squares', add(mul(c, c), mul(c2, c2)), False)
    # negative controls: random core polynomials (E1 should fail)
    for d in (6, 7, 8):
        for par in (0, 1):
            D = {m: rng.randrange(1, p) for m in W.mons(d, W.CORE, par) if rng.random() < 0.5}
            addT(f'random core-only degree {d} par{par} (control, expect E1 fail)', 'random-control', D, True)
    log(JOB, {'exp': 'struct-plan', 'n_tasks': len(tasks), 'workers': NW})
    with Pool(NW) as pool:
        for res in pool.imap_unordered(task, tasks):
            log(JOB, res)
    log(JOB, {'exp': 'done', 'seconds': round(time.time() - t0, 1)})
