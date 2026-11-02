"""Job B (excl): degree-5 denominators.  By the E1-kernel computation every E1-admissible D of degree 5 lies in
J(F)_5 = {F_v = sum_a v_a dF/dz_a : v in C^6}.  For D = F_v the general solution of E1 is
    N = A v + (m/6) z + K,   K = sum_{i<j} l_ij (F_j d_i - F_i d_j)  (Koszul, l_ij linear forms),
and E2 reduces (m drops out) to
    F_v A_v - A F_vv + F_v div K - K(F_v)  =  F Q,      F_vv = v^T Hess(F) v, A_v = v(A).
HESSIAN RELAXATION: replace the rank-one tensors v (x) v and v (x) l by arbitrary tensors S in Sym^2(C^6),
T in C^6 (x) C^{15x6}; the resulting LINEAR system is solved EXACTLY over Q; the projection W of its kernel to
Sym^2 must contain v (x) v for every certificate.  If W contains no nonzero rank-one symmetric tensor (a
Groebner computation in 6 variables over Q) then no certificate of degree 5 exists over C.
Also: direct fixed-D tests of the six partials and random rational members (pure and mixed parity) with the
exact-nullity [RANK] argument.  Outputs remote/<job>.jsonl, remote/<job>-relax.json"""
import sys, time, json, random, itertools, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
from w2lib import log, save_json
from hodge import core
from hodge.core import monomials, mul, add, sub, shift, pderiv, scale, NV
from fractions import Fraction
from flint import fmpz_mat, fmpq_mat

JOB = sys.argv[1] if len(sys.argv) > 1 else 'excl-d5'
F, A = W.F, W.A
Fi = [pderiv(F, i) for i in range(NV)]
Fij = [[pderiv(Fi[i], j) for j in range(NV)] for i in range(NV)]
Ai = [pderiv(A, i) for i in range(NV)]
Z = [{tuple(1 if j == i else 0 for j in range(NV)): Fraction(1)} for i in range(NV)]

def F_v(v):
    out = {}
    for a in range(NV):
        if v[a]: out = add(out, scale(Fi[a], Fraction(v[a])))
    return out

# ---------- 0. algebraic sanity check of the reduction (exact identity on random data) ----------
rng = random.Random(5)
def rand_poly(deg, par=None, rng=rng, lo=-3, hi=3):
    return {m: Fraction(rng.randint(lo, hi)) for m in monomials(deg, par) if rng.random() < 0.3}
t0 = time.time()
v = [Fraction(rng.randint(-3, 3)) for _ in range(NV)]
D = F_v(v)
l = {(i, j): rand_poly(1) for i in range(NV) for j in range(i + 1, NV)}
m = rand_poly(5)
N = [scale(A, v[i]) for i in range(NV)]
N = [add(N[i], scale(mul(m, Z[i]), Fraction(1, 6))) for i in range(NV)]
for (i, j), lij in l.items():
    N[i] = add(N[i], mul(lij, Fi[j])); N[j] = sub(N[j], mul(lij, Fi[i]))
E = core.defect(F, A, D, 1, N)
# predicted: E = F * (F_v A_v - A F_vv + F_v div K - K(F_v))
Av = {}; Fvv = {}
for a in range(NV):
    Av = add(Av, scale(Ai[a], v[a]))
    for b in range(NV):
        Fvv = add(Fvv, scale(Fij[a][b], v[a] * v[b]))
divK = {}; KFv = {}
Dv = [pderiv(D, i) for i in range(NV)]
for (i, j), lij in l.items():
    divK = add(divK, sub(mul(pderiv(lij, i), Fi[j]), mul(pderiv(lij, j), Fi[i])))
    KFv = add(KFv, mul(lij, sub(mul(Fi[j], Dv[i]), mul(Fi[i], Dv[j]))))
pred = mul(F, add(sub(mul(D, Av), mul(A, Fvv)), sub(mul(D, divK), KFv)))
ok_identity = (sub(E, pred) == {})
log(JOB, {'exp': 'd5-reduction-identity-check', 'result': 'identity holds' if ok_identity else 'IDENTITY FAILS',
          'seconds': round(time.time() - t0, 1), 'note': 'E(A,D=F_v,N=Av+mz/6+K) == F*(F_v A_v - A F_vv + F_v divK - K(F_v)) on random exact data'})
assert ok_identity
if os.environ.get("W2_STOP"): sys.exit(0)

# ---------- 1. the relaxed linear system, exact over Q ----------
t0 = time.time()
cols, info = [], []
# S-columns (symmetric tensor coordinates s_ab, a<=b)
for a in range(NV):
    for b in range(a, NV):
        if a == b:
            col = sub(mul(Fi[a], Ai[a]), mul(A, Fij[a][a]))
        else:
            col = sub(add(mul(Fi[a], Ai[b]), mul(Fi[b], Ai[a])), scale(mul(A, Fij[a][b]), 2))
        cols.append(col); info.append(('S', a, b))
nS = len(cols)
# T-columns: T_{a,(ij),c} ~ v_a * l_{ij,c}
for a in range(NV):
    for i in range(NV):
        for j in range(i + 1, NV):
            for c in range(NV):
                col = {}
                if c == i: col = add(col, mul(Fi[a], Fi[j]))
                if c == j: col = sub(col, mul(Fi[a], Fi[i]))
                col = sub(col, mul(Z[c], sub(mul(Fi[j], Fij[a][i]), mul(Fi[i], Fij[a][j]))))
                cols.append(col); info.append(('T', a, i, j, c))
nT = len(cols) - nS
# Q-columns: -F m
for mq in monomials(4):
    cols.append(scale(shift(F, mq, 1), -1)); info.append(('Q', mq))
nQ = len(cols) - nS - nT
M, rows, scales = W.build_exact(cols)
X, nul = M.nullspace()
rank = M.ncols() - nul
# project kernel to S-coordinates (undo column scaling: coordinate j of the true kernel vector = X[j,k]/scales[j]...
# scaling columns by s_j means M' = M diag(s); kernel of M' = diag(s)^{-1} kernel(M); so true S-coords = X[j,k]/s_j)
Wrows = []
for k in range(nul):
    vec = [Fraction(int(X[j, k]), scales[j]) for j in range(nS)]
    if any(vec): Wrows.append(vec)
if Wrows:
    den = 1
    for vec in Wrows:
        for x in vec:
            den = den * x.denominator // __import__('math').gcd(den, x.denominator)
    Wm = fmpz_mat(len(Wrows), nS, [int(x * den) for vec in Wrows for x in vec])
    R, dd, rW = Wm.rref()
    Wbasis = [[Fraction(int(R[i, j]), int(dd)) for j in range(nS)] for i in range(rW)]
else:
    rW = 0; Wbasis = []
log(JOB, {'exp': 'd5-hessian-relaxation', 'prime': None, 'rows': M.nrows(), 'cols': M.ncols(), 'nS': nS, 'nT': nT, 'nQ': nQ,
          'rank': rank, 'nullity': nul, 'dim_W': rW, 'seconds': round(time.time() - t0, 1),
          'note': 'EXACT over Q: kernel of [Psi(S) + Phi(T) - F Q] and its projection W to Sym^2(C^6)'})

# ---------- 2. rank-one tensors in W ----------
t0 = time.time()
Sidx = [(a, b) for a in range(NV) for b in range(a, NV)]
result = {}
if rW == 0:
    result = {'rank_one_points': 'none (W = 0)', 'quadrics': 0}
    log(JOB, {'exp': 'd5-rank-one-test', 'result': 'W = 0: no certificate of degree 5 exists over C [EXACT]', 'seconds': 0.0})
else:
    # annihilator of W: solve Wbasis * lambda = 0  (lambda in Q^21)
    den = 1
    for vec in Wbasis:
        for x in vec: den = den * x.denominator // __import__('math').gcd(den, x.denominator)
    Wm = fmpz_mat(rW, nS, [int(x * den) for vec in Wbasis for x in vec])
    Lam, nlam = Wm.nullspace()
    quadrics = []
    vars_ = [f'v{a}' for a in range(NV)]
    for k in range(nlam):
        terms = []
        for jj, (a, b) in enumerate(Sidx):
            c = int(Lam[jj, k])
            if c: terms.append(f'({c})*v{a}*v{b}')
        quadrics.append('+'.join(terms) if terms else '0')
    # Singular: dimension of the affine cone defined by the quadrics (0 => only v = 0)
    sing = 'ring r=0,(' + ','.join(vars_) + '),dp; ideal I=' + ','.join(quadrics) + '; ideal G=std(I); print(dim(G)); print(size(G)); print(G); quit;'
    with open(f'{W.REMOTE_OUT}/{JOB}-rankone.sing', 'w') as f: f.write(sing)
    try:
        out = subprocess.run(['Singular', '-q', f'{W.REMOTE_OUT}/{JOB}-rankone.sing'], capture_output=True, text=True, timeout=3600).stdout
        dim_cone = int(out.strip().splitlines()[0])
    except Exception as ex:
        out = f'Singular failed: {ex}'; dim_cone = None
    result = {'dim_W': rW, 'n_quadrics': nlam, 'affine_cone_dim': dim_cone, 'singular_output_head': out[:2000]}
    verdict = ('no nonzero rank-one tensor in W: no certificate of degree 5 exists over C [EXACT]' if dim_cone == 0
               else f'W contains rank-one tensors (cone dim {dim_cone}): candidates must be tested')
    log(JOB, {'exp': 'd5-rank-one-test', 'result': verdict, 'dim_W': rW, 'n_quadrics': nlam, 'cone_dim': dim_cone,
              'seconds': round(time.time() - t0, 1)})
save_json(JOB, 'relax', {'nS': nS, 'nT': nT, 'nQ': nQ, 'rank': rank, 'nullity': nul, 'dim_W': rW,
                         'W_basis': [[str(x) for x in vec] for vec in Wbasis], 'S_index': Sidx, 'result': result})

# ---------- 3. direct fixed-D tests (partials, random rational members), [RANK] ----------
names = ['dF/dX0', 'dF/dX1', 'dF/dX2', 'dF/dY0', 'dF/dY1', 'dF/dY2']
for p in (W.P2, W.P3):
    for i in range(NV):
        r = W.fixed_D_test(Fi[i], p)
        log(JOB, {'exp': 'd5-partial', 'D': names[i], 'd': 5, 'prime': p, **r})
rng = random.Random(11)
for trial in range(3):
    ve = [rng.randint(-5, 5) for _ in range(3)] + [0, 0, 0]
    vo = [0, 0, 0] + [rng.randint(-5, 5) for _ in range(3)]
    vm = [rng.randint(-5, 5) for _ in range(6)]
    for tag, v, pf in (('Y-even v', ve, False), ('Y-odd v', vo, False), ('mixed v', vm, True)):
        D = F_v(v)
        if not D: continue
        for p in (W.P2, W.P3):
            r = W.fixed_D_test(D, p, parity_free=pf)
            log(JOB, {'exp': 'd5-random-Jacobian-member', 'D': f'F_v, v={v}', 'family': tag, 'd': 5, 'prime': p, **r})
log(JOB, {'exp': 'done'})
