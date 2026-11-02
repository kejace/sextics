"""JOB excl-dp-ann: annihilator catalogue for the sum principle, and Hessian-annihilators of Jacobian-type D.

For a homogeneous D (any degree) the served space satisfies
    served(D) / J_6  subset  Ann_6(D) := { A' in R_6^- : [A'][D] = 0 in R }        (E1 for A')
and by the sum principle served(D1 D2) = served(D1) + served(D2) (classes supported on V(D1) u V(D2)).
Hence  A in served(prod D_i)  ==>  [A] in sum_i Ann_6(D_i)  (mod J).   Kernel dimensions mod p bound the rational
ones from above (modular rank <= rational rank), so  A notin sum Ann_6^p(D_i)  is a proof over Q.
For D = w(F) in J (vacuous Ann), the relaxation gives served(D)/J_6 subset Ann_6([H(w,w)]).

Candidates: all structured low-degree factors used in the dossier (§4.2) and wave 1, quadrics/cubics/quartics with
symmetry, orbit sums O1,O2,O3, Vandermonde squares, coupling part, plus Hessian forms of the six partials and the
symmetric fields.  Output: ann_<job>.json; log lines exp='ann'."""
import sys, os, time, json
import numpy as np
import w2lib as W
from hodge import core
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree, parse

job = os.environ.get('W2JOB', 'ann')
OUT = f'{W.PKT}/ann_{job}.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
p = W.P1
F, A = W.F, W.A
Ap = reduce_mod(A, p)

X = ['X0', 'X1', 'X2']; Y = ['Y0', 'Y1', 'Y2']
sym = lambda s: s
cands = {
    # degree 2
    'p2X': 'X0**2+X1**2+X2**2', 'p2Y': 'Y0**2+Y1**2+Y2**2', 'pXY': 'X0*Y0+X1*Y1+X2*Y2',
    'e2X': 'X0*X1+X1*X2+X2*X0', 'e2Y': 'Y0*Y1+Y1*Y2+Y2*Y0',
    'X0^2+Y0^2': 'X0**2+Y0**2', 'X0^2-Y0^2': 'X0**2-Y0**2', 'X0X1+Y0Y1': 'X0*X1+Y0*Y1', 'X0X1-Y0Y1': 'X0*X1-Y0*Y1',
    '(X0-X1)^2+(Y0-Y1)^2': '(X0-X1)**2+(Y0-Y1)**2', '(X0-X1)^2-(Y0-Y1)^2': '(X0-X1)**2-(Y0-Y1)**2',
    'X0Y1-X1Y0': 'X0*Y1-X1*Y0', 'X0Y1+X1Y0': 'X0*Y1+X1*Y0', 'X2^2': 'X2**2', 'X0X1': 'X0*X1', 'X2Y2': 'X2*Y2',
    # degree 3
    'p3X': 'X0**3+X1**3+X2**3', 'e3X': 'X0*X1*X2', 'VdmX': '(X0-X1)*(X1-X2)*(X2-X0)', 'VdmY': '(Y0-Y1)*(Y1-Y2)*(Y2-Y0)',
    'sXY2': 'X0*Y0**2+X1*Y1**2+X2*Y2**2', 'sX2Y': 'X0**2*Y0+X1**2*Y1+X2**2*Y2', 'XYYcyc': 'X0*Y1*Y2+X1*Y0*Y2+X2*Y0*Y1',
    'Hesse': 'X0**3+X1**3+X2**3-3*X0*X1*X2', 'p3Y': 'Y0**3+Y1**3+Y2**3', 'e3Y': 'Y0*Y1*Y2',
    # degree 4
    'p2X^2': '(X0**2+X1**2+X2**2)**2', 'p2Xp2Y': '(X0**2+X1**2+X2**2)*(Y0**2+Y1**2+Y2**2)', 'pXY^2': '(X0*Y0+X1*Y1+X2*Y2)**2',
    'p4X': 'X0**4+X1**4+X2**4', 'sX2Y2': 'X0**2*Y0**2+X1**2*Y1**2+X2**2*Y2**2', 'p4Y': 'Y0**4+Y1**4+Y2**4',
    'X0^4-X0^2Y0^2+Y0^4': 'X0**4-X0**2*Y0**2+Y0**4', 'X0^4-X1^2Y0^2': 'X0**4-X1**2*Y0**2',
    # degree 6 (non-Jacobian ones; the orbit sums and Vandermonde squares)
    'VdmX^2': '((X0-X1)*(X1-X2)*(X2-X0))**2', 'VdmY^2': '((Y0-Y1)*(Y1-Y2)*(Y2-Y0))**2',
    'VdmXVdmY': '(X0-X1)*(X1-X2)*(X2-X0)*(Y0-Y1)*(Y1-Y2)*(Y2-Y0)',
    'prod+': '((X0-X1)**2+(Y0-Y1)**2)*((X1-X2)**2+(Y1-Y2)**2)*((X2-X0)**2+(Y2-Y0)**2)',
    'prod-': '((X0-X1)**2-(Y0-Y1)**2)*((X1-X2)**2-(Y1-Y2)**2)*((X2-X0)**2-(Y2-Y0)**2)',
    'p6X': 'X0**6+X1**6+X2**6', 'p6Y': 'Y0**6+Y1**6+Y2**6', 'sumT': 'X0**3*Y0**3+X1**3*Y1**3+X2**3*Y2**3',
    'B2': 'X0**3*Y0**3+X1**3*Y1**3-2*X2**3*Y2**3', 'A': 'X0**3*Y0**3-X1**3*Y1**3', 'pXY^3': '(X0*Y0+X1*Y1+X2*Y2)**3',
}
# orbit sums O1, O2, O3 and the coupling part
O1 = []; O2 = []; O3 = []
for (i, j, k) in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
    O1.append(f'X{i}**3*X{j}*Y{i}*Y{j}'); O2.append(f'X{i}**3*X{j}*Y{k}**2'); O3.append(f'Y{i}**3*Y{j}*X{k}**2')
cands['O1'] = '+'.join(O1); cands['O2'] = '+'.join(O2); cands['O3'] = '+'.join(O3)
cands['Fc'] = f'5*({cands["O1"]})/7+5*({cands["O2"]})/11+5*({cands["O3"]})/13'


def save():
    json.dump(res, open(OUT, 'w'), indent=1, default=str)


# reducers for R_6^- and the targets R_{6+d}
red6 = W.Red(6, p, 1, None, verbose=True)
std6 = red6.std_mons
reds = {}


def ann6(Dp, label, extra=None):
    """Ann_6(D) = {A' in R_6^- : A' D in J}; returns (dim, contains A, basis as polys)."""
    d = degree(Dp)
    par = core.parity_of(Dp)
    assert par is not None
    key = (6 + d, (1 + par) % 2)
    if key not in reds:
        reds[key] = W.Red(key[0], p, key[1], None, verbose=True)
    redT = reds[key]
    t0 = time.time()
    prods = [mul({m: 1}, Dp, p) for m in std6]
    M = redT.nf_polys(prods)
    ker = W.nullspace_np(np.ascontiguousarray(M.T), p)
    dim = ker.shape[0]
    # does A lie in the kernel (mod J)?  A's coordinates in std6:
    Av = red6.nf(Ap)
    contains = False
    if dim:
        st = np.vstack([ker, Av.reshape(1, -1)])
        contains = (W.rank_np(st, p) == dim)
    out = {'label': label, 'd': d, 'par': par, 'dim_Ann6': int(dim), 'contains_A': bool(contains),
           'rank': int(M.shape[0] - dim), 'seconds': round(time.time() - t0, 1)}
    if extra:
        out.update(extra)
    return out, ker


mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
if mode in ('all', 'cands'):
    for label, expr in cands.items():
        if label in res:
            continue
        D = parse(expr)
        Dp = reduce_mod(D, p)
        inJ = (degree(D) >= 5) and W.jacobian.in_jacobian(F, D, p)[0]
        if inJ:
            res[label] = {'label': label, 'd': degree(D), 'in_J': True, 'note': 'Jacobian-type: Ann vacuous'}
            save(); continue
        out, ker = ann6(Dp, label, {'expr': expr, 'in_J': False})
        res[label] = out; save()
        W.log('ann', D=label, d=out['d'], prime=p, result=f"dim Ann6={out['dim_Ann6']} containsA={out['contains_A']}",
              rank=out['rank'], seconds=out['seconds'], note='Ann_6(D) = {A\' in R_6^- : A\' D in J}; served(D)/J_6 inside it')
        print(label, out, flush=True)

if mode in ('all', 'hess'):
    # Hessian forms of structured fields w (D = w(F) in J): served(w(F))/J_6 subset Ann_6([H(w,w)])
    Fp, Fi, Fij = W.hessian(p)
    fields = {}
    for i in range(6):
        w = [{} for _ in range(6)]; w[i] = {(0,) * 6: 1}; fields[f'd_{i}'] = w
    def lin(i, j):
        e = [0] * 6; e[j] = 1; return {tuple(e): 1}
    w = [{} for _ in range(6)]
    for i in range(3): w[i] = lin(i, i)
    fields['EulerX'] = w
    w = [{} for _ in range(6)]; w[0] = lin(0, 0); w[3] = lin(3, 3); fields['pairradial0'] = w
    w = [{} for _ in range(6)]; w[0] = lin(0, 0); w[3] = lin(3, 3); w[1] = scale(lin(1, 1), -1, p); w[4] = scale(lin(4, 4), -1, p); fields['gA'] = w
    w = [{} for _ in range(6)]; w[0] = lin(0, 1); w[1] = lin(1, 0); w[3] = lin(3, 4); w[4] = lin(4, 3); fields['swap01'] = w
    w = [{} for _ in range(6)]; w[0] = lin(0, 1); w[1] = scale(lin(1, 0), -1, p); w[3] = lin(3, 4); w[4] = scale(lin(4, 3), -1, p); fields['rot01'] = w
    w = [{} for _ in range(6)]
    for i in range(3): w[i] = lin(i, (i + 1) % 3); w[3 + i] = lin(3 + i, 3 + (i + 1) % 3)
    fields['cycshift'] = w
    w = [{} for _ in range(6)]; w[0] = {(0,) * 6: 1}; w[1] = {(0,) * 6: 1}; w[2] = {(0,) * 6: 1}; fields['d0+d1+d2'] = w
    w = [{} for _ in range(6)]; w[0] = {(0,) * 6: 1}; w[1] = {(0,) * 6: p - 1}; fields['d0-d1'] = w
    for label, w in fields.items():
        key = f'hess_{label}'
        if key in res:
            continue
        Hw = {}
        for i in range(6):
            for j in range(6):
                if w[i] and w[j]:
                    Hw = add(Hw, mul(mul(w[i], w[j], p), Fij[i][j], p), p)
        out, ker = ann6(Hw, key, {'field': label, 'deg_D': 5 + degree(w[next(i for i in range(6) if w[i])])})
        res[key] = out; save()
        W.log('ann_hess', D=f'w(F) with w={label}', d=out['deg_D'], prime=p,
              result=f"dim Ann6(Hw)={out['dim_Ann6']} containsA={out['contains_A']}", rank=out['rank'], seconds=out['seconds'],
              note='relaxation: served(w(F))/J_6 inside Ann_6([H(w,w)])')
        print(key, out, flush=True)
print('DONE', flush=True)
