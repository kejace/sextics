"""Job: structured candidate sweep (exact sector solver; E1 screen first). usage: python3 w2_sweep.py [procs]"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import time, json, itertools, random
import multiprocessing as mp
from fractions import Fraction
from w2core import *
from hodge import jacobian
procs = int(sys.argv[1]) if len(sys.argv) > 1 else 8
JOB = 'w2sweep'
F, A, q = load()
P = core.parse
M = core.mul
t = Fraction(5, 7)
V = core.VARS
X0, X1, X2, Y0, Y1, Y2 = [P(v) for v in V]
C = P('X0**3 + X0*Y1*Y2 + X1*Y0*Y1')
H = P('X0**6 + X1**6 + Y0**6 + Y1**6 + Y2**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0**2*Y0*Y1*Y2**2/7 + 5*X0*X1*Y1**3*Y2/7')
Z = P('X2**3'); Zp = core.add(core.scale(Z, 2), core.scale(C, t)); Hp = core.sub(H, core.scale(M(C, C), t * t / 4))
Fi = [pderiv(F, i) for i in range(6)]
D8 = P(json.load(open(f'{WD}/kernel_d8_0_4_0.json'))['D'])
cands = []
def addc(D, name, fam, note=''):
    if D and core.is_homogeneous(D) and core.parity_of(D) is not None:
        cands.append((D, name, fam, note))
def restrict(Pl, zero):
    return {e: v for e, v in Pl.items() if all(e[i] == 0 for i in zero)}
# (1) restrictions of F to coordinate subspaces: cones with vertex a point / line / plane
for k in (1, 2, 3):
    for zero in itertools.combinations(range(6), k):
        D = restrict(F, zero)
        nm = 'F|_{' + ','.join(V[i] for i in zero) + '=0}'
        addc(D, nm, f'F-restriction-{k}', f'cone with vertex the {["point","line","plane"][k-1]} {{{",".join(V[i] for i in zero)}}}')
# complements: the terms of F involving a given variable (F - F|_{Z_i=0}) and the coupling part
addc(core.sub(F, restrict(F, (2,))), 'F - F|_{X2=0}', 'F-pieces', 'X2^6 + X2^3 t C')
Fer = P('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
addc(Fer, 'Fermat part sum Z_i^6', 'F-pieces'); addc(core.sub(F, Fer), 'coupling part (F - Fermat)', 'F-pieces')
FX = restrict(F, (3, 4, 5)); FY = restrict(F, (0, 1, 2))
for lam in (Fraction(1), Fraction(-1), Fraction(1, 2), Fraction(2)):
    addc(core.add(FX, core.scale(FY, lam)), f'F_X + {lam} F_Y', 'join-pencil', 'deformation of the Fermat join divisor')
for lam in (Fraction(0), Fraction(1), Fraction(-5, 7), Fraction(10, 7), Fraction(1, 7)):
    addc(core.add(P('X0**6 + X1**6 + X2**6'), core.scale(P('X0**3*X2**3'), lam)), f'X0^6+X1^6+X2^6 + {lam} X0^3X2^3', 'FX-pencil')
# (2) block radial derivatives and their products (degree 12)
EX = {}; EY = {}
for i in range(3):
    EX = core.add(EX, M(P(V[i]), Fi[i])); EY = core.add(EY, M(P(V[i + 3]), Fi[i + 3]))
addc(M(EX, EY), 'E_X(F) * E_Y(F)', 'block-radial-products'); addc(M(EX, EX), 'E_X(F)^2', 'block-radial-products')
addc(M(EX, FY), 'E_X(F) * F_Y', 'block-radial-products'); addc(M(FX, EY), 'F_X * E_Y(F)', 'block-radial-products')
addc(M(FX, FY), 'F_X * F_Y', 'join-products', 'cone over Gamma_X times cone over Gamma_Y'); addc(M(FX, FX), 'F_X^2', 'join-products')
addc(M(FX, H), 'F_X * H', 'join-products'); addc(M(FY, H), 'F_Y * H', 'join-products')
# (3) Hessians of ternary/binary restrictions (degree 12 / 8 / ...)
def hess_det(Pl, vars_):
    n = len(vars_)
    import sympy as spy
    xs = spy.symbols(' '.join(V))
    expr = sum(spy.Rational(Fraction(v).numerator, Fraction(v).denominator) * spy.prod([xs[i] ** e[i] for i in range(6)]) for e, v in Pl.items())
    Hm = spy.Matrix(n, n, lambda a, b: spy.diff(expr, xs[vars_[a]], xs[vars_[b]]))
    det = spy.expand(Hm.det())
    return P(str(det).replace('**', '**'))
try:
    for vars_ in [(0, 1, 2), (3, 4, 5), (0, 2, 4), (1, 3, 4), (0, 3, 5)]:
        zero = tuple(i for i in range(6) if i not in vars_)
        D = restrict(F, zero)
        Hd = hess_det(D, vars_)
        addc(Hd, 'Hess(F|_{' + ','.join(V[i] for i in vars_) + '})', 'hessian-restriction', 'Hessian of a ternary sextic restriction (degree 12)')
except Exception as ex:
    print('hessian construction failed', ex, flush=True)
# (4) batchB2 list (never launched in wave 1)
random.seed(5)
def polar(a):
    D = {}
    for ai, Fi_ in zip(a, Fi):
        if ai: D = core.add(D, core.scale(Fi_, Fraction(ai)))
    return D
addc(polar([random.randint(-5, 5) for _ in range(6)]), 'polar F_a random a', 'polar-pencil')
for a, nm in [((0,1,0,0,0,1), 'a on line {X0=Y1=Y0=0}'), ((0,0,0,1,0,1), 'a on line {X0=Y1=X1=0}'), ((0,0,0,0,1,0), 'a = e_Y1'),
              ((1,0,0,1,0,0), 'a = e_X0+e_Y0'), ((1,0,0,-1,0,0), 'a = e_X0-e_Y0'), ((1,1,0,0,0,0), 'a=e_X0+e_X1'), ((0,0,1,1,0,0), 'a=e_X2+e_Y0')]:
    addc(polar(a), f'polar F_a, {nm}', 'polar-pencil')
for i, j in [(0, 3), (2, 2), (0, 0), (3, 3), (1, 4), (2, 5)]:
    addc(M(Fi[i], Fi[j]), f'F_{V[i]}*F_{V[j]}', 'polar-products')
for D, name, fam in [
 (M(M(X2, X2), H), 'X2^2*H', 'pieces-8'), (M(M(X2, X2), Hp), "X2^2*H'", 'pieces-8'), (M(M(X2, X2), M(C, Zp)), "X2^2*C*Z'", 'pieces-8'),
 (M(C, M(C, M(X2, X2))), 'X2^2*C^2', 'pieces-8'), (M(Zp, M(Zp, M(X2, X2))), "X2^2*Z'^2", 'pieces-8'),
 (M(X2, D8), 'X2*D8', 'D8-products'), (M(M(X2, X2), D8), 'X2^2*D8', 'D8-products'), (M(C, D8), 'C*D8', 'D8-products'), (M(Zp, D8), "Z'*D8", 'D8-products'),
 (M(Z, H), 'X2^3*H', 'pieces-9'), (M(Z, Hp), "X2^3*H'", 'pieces-9'), (M(C, Hp), "C*H'", 'pieces-9'), (M(Zp, Hp), "Z'*H'", 'pieces-9'), (M(C, M(C, C)), 'C^3', 'pieces-9'),
 (M(Zp, M(Zp, Zp)), "Z'^3", 'pieces-9'), (M(C, M(C, Zp)), "C^2*Z'", 'pieces-9'), (M(C, M(Zp, Zp)), "C*Z'^2", 'pieces-9'),
 (M(H, M(C, Zp)), "H*C*Z'", 'pieces-12'), (M(H, Hp), "H*H'", 'pieces-12'), (M(H, M(C, C)), 'H*C^2', 'pieces-12'), (M(Hp, M(C, Zp)), "H'*C*Z'", 'pieces-12'),
 (M(H, M(Zp, Zp)), "H*Z'^2", 'pieces-12'), (M(H, H), 'H^2', 'pieces-12'), (M(Hp, Hp), "H'^2", 'pieces-12'), (M(D8, D8), 'D8^2', 'D8-products'),
 (M(FX, D8), 'F_X*D8', 'D8-products'), (M(FY, D8), 'F_Y*D8', 'D8-products')]:
    addc(D, name, fam)
print(f'{len(cands)} candidates', flush=True)

def run(args):
    D, name, fam, note = args
    t0 = time.time()
    d = core.degree(D)
    rec = {'batch': 'sweep', 'job': JOB, 'family': fam, 'name': name, 'd': d, 'terms': len(D), 'note': note}
    try:
        try:
            divF, _ = jacobian.divisible_by(F, D) if d >= 6 else (False, None)
        except Exception:
            divF = None
        rec['F_divides_D'] = divF
        if divF:
            rec['consistent'] = None; rec['skip'] = 'F | D'; rec['seconds'] = round(time.time() - t0, 1); return rec
        ok, info = fixed_d.e1_solvable(F, A, q, D)
        rec['e1'] = bool(ok); rec['e1_info'] = {k: info[k] for k in ('rank', 'rank_aug', 'rows', 'cols')}
        if not ok:
            rec['consistent'] = False; rec['seconds'] = round(time.time() - t0, 1); return rec
        res = sect.test(F, A, q, D, exact=(d <= 9), verbose=False, lift=True)
        rec.update({k: v for k, v in res.items() if k not in ('N', 'Q')})
        if res.get('verified'):
            rec['WITNESS'] = True
            with open(f'{WD}/witness-candidate.json', 'w') as f:
                json.dump({'instanceId': IID, 'job': JOB, 'name': name, 'd': d, 'D': core.to_str(D), 'N': res['N'], 'Q': res['Q']}, f, indent=1)
    except Exception as ex:
        rec['error'] = repr(ex)
    rec['seconds'] = round(time.time() - t0, 1)
    rec['D'] = core.to_str(D)[:300]
    return rec

cands.sort(key=lambda c: core.degree(c[0]))
with mp.get_context('fork').Pool(procs) as pool:
    for rec in pool.imap_unordered(run, cands):
        w2log(rec)
        print(f"d={rec['d']:2d} e1={rec.get('e1')} cons={rec.get('consistent')} rank={rec.get('rank')}/{rec.get('rank_aug')} "
              f"{rec.get('seconds')}s {rec['family']:22s} {rec['name'][:60]} {'*** WITNESS ***' if rec.get('WITNESS') else ''} {rec.get('error','')}", flush=True)
print('done', flush=True)
