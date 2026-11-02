"""Structured irreducible candidate denominators of degree 8..12 for the E1 screen (tag<TAB>expr)."""
import sys, itertools, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import mul, add, sub, scale, power, pderiv, substitute_scalar
from fractions import Fraction
F, A, q, rec = core.load_instance('reflection-coupled-residual')
V = core.VARS
out = []
def emit(tag, P):
    if P and core.is_homogeneous(P) and core.parity_of(P) is not None:
        out.append((tag, core.to_str(P)))
# discriminant of F as a depressed cubic in u = Y2^2:  F = u^3 + Pu*u + Qc
Qc = substitute_scalar(F, 5, 0)
Pu = {}
for e, c in F.items():
    if e[5] == 2:
        e2 = list(e); e2[5] = 0; Pu[tuple(e2)] = c
Delta = sub(scale(power(Pu, 3), -4), scale(power(Qc, 2), 27))
emit('disc-u-d12', Delta)
emit('disc-u-part-P3-d12', power(Pu, 3))
emit('Qc2-d12', power(Qc, 2))   # reducible control
# branch-type combos: Qc * u-part etc. are reducible; skip.  Resultant-type: Pu^3 + t Qc^2 one-parameter family
for t in (1, -1, 2, -2, 4, -4, 27, -27):
    emit(f'P3+{t}Qc2-d12', add(power(Pu, 3), scale(power(Qc, 2), t)))
# Fermat-type cones (vertex a plane / a line) at degrees 8,10,12 and pure-Y odd degrees 9,11
def fermat(vars_, k, coeffs=None):
    P = {}
    for j, v in enumerate(vars_):
        e = [0]*6; e[v] = k
        P[tuple(e)] = Fraction(coeffs[j] if coeffs else 1)
    return P
for k in (8, 10, 12):
    for tri in itertools.combinations(range(6), 3):
        emit(f'fermat3-{"".join(V[i] for i in tri)}-d{k}', fermat(tri, k))
    for quad in itertools.combinations(range(6), 4):
        emit(f'fermat4-{"".join(V[i] for i in quad)}-d{k}', fermat(quad, k))
for k in (9, 11):
    for tri in itertools.combinations((3, 4, 5), 3):
        emit(f'fermatY3-d{k}', fermat(tri, k))
    for tri in itertools.combinations((0, 1, 2), 3):
        emit(f'fermatX3-d{k}', fermat(tri, k))
    for pair in itertools.combinations((3, 4, 5), 2):
        emit(f'fermatY2-{"".join(V[i] for i in pair)}-d{k}', fermat(pair, k))
# mixed-sign / coupled ternary forms built from the instance's own blocks at degree 8 and 12
X03Y03 = A
emit('A-squared-control-d12', power(A, 2))
emit('Qc-Pu2-d10?', mul(Qc, power(Pu, 1)))   # degree 10 reducible control
# accepted explicit witnesses with the same target on other varieties (transport test)
for iid in ('family-hesse-web', 'family-block-cube'):
    try:
        d, D, N, Q = core.load_witness(iid)
        emit(f'transport-{iid}-d{d}', D)
        emit(f'transport-{iid}-sq-d{2*d}', power(D, 2))
    except Exception as ex:
        print('skip', iid, ex)
with open('cands_w2.txt', 'w') as f:
    for tag, expr in out:
        f.write(f'{tag}\t{expr}\n')
print(len(out), 'candidates written')
