"""Validate the own FixedD builder: (a) random (N,Q) vector times matrix == direct polynomial evaluation of the
integer identity; (b) accepted witnesses fermat-q2 / block-q2 are consistent, controls inconsistent."""
import sys, json, random, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify')
import vcheck as V
from fractions import Fraction
F, A, q = V.load_instance()
# (a) matrix self-test on cycle-q2 with D = F_{X2}
D = V.pderiv(F, 2); D = V.scale(D, 2)            # 2*F_{X2} integral
S = V.FixedD(F, A, q, D, verbose=True)
random.seed(7)
x = [random.randint(-5, 5) for _ in range(S.ncols)]
N = [dict() for _ in range(6)]; Q = {}
for (kind, i, m), v in zip(S.meta, x):
    if v == 0: continue
    if kind == 'N': N[i][m] = Fraction(v)
    else: Q[m] = Fraction(v)
F2 = V.scale(F, 2); Dp = {m: Fraction(v) for m, v in D.items()}
c = 2 ** q
E = V.scale(V.mul(A, V.mul(Dp, Dp)), 2 ** (q + 1))
for i in range(6):
    E = V.add(E, V.mul(Dp, V.mul(N[i], V.pderiv(F2, i))), -c * q)
    E = V.add(E, V.mul(F2, V.add(V.mul(Dp, V.pderiv(N[i], i)), V.mul(N[i], V.pderiv(Dp, i)), -1)), c)
E = V.add(E, V.mul(V.power(F2, q + 1), Q), -1)
# matrix: M x - b should equal E  (system is M x = b with b = -2^(q+1) A D^2 ... check sign: cols sum + rhs_col = E)
Mx = {}
for j, col in enumerate(S.cols):
    if x[j]: Mx = V.add(Mx, col, x[j])
E_mat = V.add(Mx, S.rhs, -1)      # M x - b, where S.rhs = -2^(q+1) A D^2 => M x + 2^(q+1) A D^2
ok = (E_mat == {m: Fraction(v) for m, v in E.items()} or {m: Fraction(v) for m, v in E_mat.items()} == {m: Fraction(v) for m, v in E.items()})
print('matrix self-test ok:', ok, 'terms', len(E))
V.log(claim='own FixedD matrix equals direct evaluation of the integer identity (random N,Q; D=2F_X2)', source='self', verdict='CONFIRMED' if ok else 'REFUTED', evidence={'terms': len(E), 'cols': S.ncols, 'rows': S.nrows})
# (b) accepted witnesses
for iid in ['fermat-q2', 'block-q2']:
    Fi, Ai, qi = V.load_instance(f'/home/kejace/devel/sextics/data/{iid}/instance.rec')
    w = json.load(open(f'/home/kejace/devel/sextics/data/{iid}/witness.json'))
    Di = V.parse(w['D'])
    # make integral
    import math
    den = 1
    for v in Di.values(): den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    Di = V.scale(Di, den)
    t0 = time.time()
    Si = V.FixedD(Fi, Ai, qi, Di, verbose=True)
    r0, ra = Si.ranks_mod(999983); e0, ea = Si.ranks_exact()
    Dc = V.add(Di, {(V.degree(Di), 0, 0, 0, 0, 0): Fraction(1)})
    Sc = V.FixedD(Fi, Ai, qi, Dc)
    c0, ca = Sc.ranks_mod(999983)
    res = {'instance': iid, 'D': V.to_str(Di), 'mod_999983': (r0, ra), 'exact': (e0, ea), 'control_mod': (c0, ca), 'seconds': round(time.time() - t0, 1)}
    print(res)
    V.log(claim=f'own FixedD solver: accepted witness {iid} consistent, control D+X0^d inconsistent', source='self-validation', verdict='CONFIRMED' if (r0 == ra and e0 == ea and c0 + 1 == ca) else 'REFUTED', evidence=res)
