"""Second-prime / second-method re-verification of wave-1 modular claims.
usage: python3 w2x_reverify_job.py --prime P --parts engine,twisted,jinside,invker,served"""
import sys, time, json, argparse, random
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
import w2x_lib as L
import charsolve as CS
from hodge import core, fixed_d, jacobian
from fractions import Fraction
ap = argparse.ArgumentParser()
ap.add_argument('--prime', type=int, default=L.P2); ap.add_argument('--parts', default='twisted,engine,jinside,invker,served')
a = ap.parse_args()
F, A, q = L.load(); p = a.prime
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
X, Y = L.X, L.Y
DX, DY = L.dlt(X), L.dlt(Y)
e = {('X', k): L.e_k(X, k) for k in (1, 2, 3)}; e.update({('Y', k): L.e_k(Y, k) for k in (1, 2, 3)})
QX, QY, e1X, e1Y, e3X, e3Y = e[('X', 2)], e[('Y', 2)], e[('X', 1)], e[('Y', 1)], e[('X', 3)], e[('Y', 3)]
e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e1Y, e3Y), 2))
D7, D7Y = mul(DX, e2sqX), mul(DY, e2sqY)
Fi = [core.pderiv(F, i) for i in range(6)]
def radial(idx):
    out = {}
    for i in idx: out = add(out, mul(L.var(i), Fi[i]))
    return out
EXF, EYF = radial((0, 1, 2)), radial((3, 4, 5))
GG, GsX, GsY, Gss = L.group('triv'), L.group('sgnX'), L.group('sgnY'), L.group('sgnsgn')
parts = a.parts.split(',')
def twisted(label, D, G, fam):
    t0 = time.time()
    S = CS.CharSystem(F, A, q, D, G, p=p); r = S.solve()
    rec = {'exp': 'reverify-twisted', 'D': label, 'd': core.degree(D), 'prime': p, 'family': fam, 'result': 'consistent' if r['consistent'] else 'inconsistent',
           'rank': r['rank'], 'nullity': r['nullity'], 'rows': S.nrows, 'cols': S.ncols, 'seconds': round(time.time() - t0, 1), 'job': 'reverify'}
    L.log(rec); print(json.dumps(rec), flush=True); return rec
def engine(label, D, fam):
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
    if not ok:
        rec = {'exp': 'reverify-engine', 'D': label, 'd': core.degree(D), 'prime': p, 'family': fam, 'result': 'E1 fails', 'rank': info.get('rank'), 'seconds': round(time.time() - t0, 1), 'job': 'reverify'}
    else:
        r = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, want_exact=True)
        rec = {'exp': 'reverify-engine', 'D': label, 'd': core.degree(D), 'prime': p, 'family': fam, 'result': 'consistent' if r['consistent'] else 'inconsistent',
               'rank': r.get('rank'), 'nullity': r.get('nullity'), 'rows': r.get('rows'), 'cols': r.get('cols'), 'exact': bool(r.get('exact')), 'seconds': round(time.time() - t0, 1), 'job': 'reverify'}
    L.log(rec); print(json.dumps(rec), flush=True); return rec
if 'twisted' in parts:
    twisted('D7=DX*e2(X^2)', D7, GsX, 'e1-kernel-deg7'); twisted('D7Y=DY*e2(Y^2)', D7Y, GsY, 'e1-kernel-deg7')
    twisted('EX(F)', EXF, GG, 'J6-radial'); twisted('EY(F)', EYF, GG, 'J6-radial')
    twisted('DX*DY*QX', mul(mul(DX, DY), QX), Gss, 'sgnsgn8'); twisted('DX*DY*e1X', mul(mul(DX, DY), e1X), Gss, 'sgnsgn7')
if 'engine' in parts:
    engine('D7=DX*e2(X^2)', D7, 'e1-kernel-deg7'); engine('D7Y=DY*e2(Y^2)', D7Y, 'e1-kernel-deg7')
if 'jinside' in parts:
    for lab, D in [('F_X0', Fi[0]), ('F_Y0', Fi[3]), ('F_X0+F_X1', add(Fi[0], Fi[1])), ('X0*F_X0', mul(L.var(0), Fi[0])), ('Y0*F_Y0', mul(L.var(3), Fi[3])),
                   ('X0*F_Y0', mul(L.var(0), Fi[3])), ('Y0*F_X0', mul(L.var(3), Fi[0])), ('X1*F_X0', mul(L.var(1), Fi[0])), ('EX(F)', EXF), ('EY(F)', EYF)]:
        engine(lab, D, 'J-inside-full')
if 'invker' in parts:
    random.seed(11)
    for d, pD, G, lab in [(8, 0, GG, 'triv'), (9, 0, GG, 'triv'), (9, 1, GG, 'triv'), (8, 1, GsX, 'sgnX'), (8, 1, GsY, 'sgnY'), (9, 0, GsX, 'sgnX'), (9, 1, GsY, 'sgnY')]:
        t0 = time.time()
        info, polys = CS.e1_kernel_char(F, A, q, d, pD, G, p=p)
        rec = {'exp': 'reverify-kernel-char', 'd': d, 'pD': pD, 'chi': lab, 'prime': p, 'info': info, 'seconds': round(time.time() - t0, 1), 'job': 'reverify'}
        L.log(rec); print(json.dumps(rec), flush=True)
        Jp = CS.J_char_basis(F, d, pD, G, p)
        for i, P in enumerate(polys):
            twisted(f'{lab}-ker-d{d}-p{pD}-rep{i}', core.reduce_mod(P, p), G, f'E1-kernel-{lab}-d{d}')
        for t in range(3):
            if not polys: break
            D = {}
            for P in polys: D = core.add(D, core.scale(P, random.randrange(1, p), p), p)
            for P in Jp: D = core.add(D, core.scale(P, random.randrange(1, p), p), p)
            twisted(f'{lab}-ker-d{d}-p{pD}-rand{t}(incl J)', core.reduce_mod(D, p), G, f'E1-kernel-{lab}-d{d}')
if 'served' in parts:
    import invsolve
    Treps = invsolve.orbit_reps(6, 1, invsolve.canon_codes)
    def vec(P): return [P.get(r, 0) % p for r in Treps]
    from flint import nmod_mat
    def rank_of(vs):
        if not vs: return 0
        return nmod_mat(len(vs), len(Treps), [v[i] for v in vs for i in range(len(Treps))], p).rank()
    Jv = [vec(core.reduce_mod(P, p)) for P in CS.J_char_basis(F, 6, 1, GG, p)]
    rJ = rank_of(Jv); Av = vec(core.reduce_mod(A, p))
    import scan_served
    span = list(Jv)
    for lab, D, fam in scan_served.cands:
        if fam == 'quot-curve6X': continue
        t0 = time.time()
        divF, _ = jacobian.divisible_by(F, D, p=p)
        tl, G_ = None, None
        for l2, G2 in [('triv', GG), ('sgnX', GsX), ('sgnY', GsY), ('sgnsgn', Gss)]:
            if G2.is_semi_invariant(D): tl, G_ = l2, G2; break
        if tl is None or core.parity_of(D) is None: continue
        S = CS.CharSystem(F, A, q, D, G_, p=p, targets=True)
        vs = [vec(P) for P in S.served_targets()]
        beyond = rank_of(Jv + vs) - rJ
        rec = {'exp': 'reverify-served', 'D': lab, 'd': core.degree(D), 'family': fam, 'chi': tl, 'prime': p, 'F_divides_D': bool(divF), 'served_total': len(vs), 'served_beyond_J': beyond, 'seconds': round(time.time() - t0, 1), 'job': 'reverify'}
        if not divF:
            span = span + vs
            rec['cum_beyond_J_legal'] = rank_of(span) - rJ; rec['A_in_legal_span'] = rank_of(span + [Av]) == rank_of(span)
        L.log(rec); print(json.dumps(rec), flush=True)
