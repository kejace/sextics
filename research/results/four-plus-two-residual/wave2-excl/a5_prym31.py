"""(a) Independent re-check of the Prym simplicity claim at p = 31 (and 37): brute-force point counts of
T_B : w^3 = b(t), b(t) = B(1,t), and E' : w^3 = g(v), g(v) = p(1,v) with B = p(X2^2, Y2^2), over F_p and F_{p^2},
computed with a hand-written F_{p^2} = F_p[s]/(s^2 - nonresidue) arithmetic (no shared code with wave 1).
From #T_B, #E' over F_p, F_{p^2}: traces a1, a2 of Frobenius on H^1(P) (P = Prym, dim 3), hence the power sums
s1 = a1, s2 = a2 of the 6 Frobenius eigenvalues; with the functional equation the L-polynomial of P has the form
1 + c1 T + c2 T^2 + c3 T^3 + p c2 T^4 + p^2 c1 T^5 + p^3 T^6 with c1 = -s1, c2 = (s1^2 - s2)/2; c3 is NOT determined by
F_p, F_{p^2} counts alone over Q (wave 1 used the Z[zeta] structure); we additionally count over F_{p^3} to fix c3
from s3 via Newton: c3 = -(s3 - c1 s2 ... ) [e3 = (s1^3 - 3 s1 s2 + 2 s3)/6, c3 = -e3].  Then factor over Q."""
import sys, json, itertools, time; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
import sympy
p = int(sys.argv[1]) if len(sys.argv) > 1 else 31
# b(t) = 1 + 15/7 t^2 + 10/7 t^4 + t^6 ; g(v) = 1 + 15/7 v + 10/7 v^2 + v^3
inv7 = pow(7, -1, p)
bco = [1, 0, 15 * inv7 % p, 0, 10 * inv7 % p, 0, 1]
gco = [1, 15 * inv7 % p, 10 * inv7 % p, 1]
def count_affine_cubic_cover(coeffs, field_elems, mul, add, powk, zero, one, q):
    """#{(t,w) in F_q^2 : w^3 = f(t)} + points at infinity (3 if deg f = 6 divisible by 3 with leading coeff 1 -> three points at infinity for
    the smooth model; for the cubic g: one point at infinity (deg 3, w^3 = v^3 + ...: 3 points? handled separately)"""
    pass
# generic finite field F_{p^k} as polynomials mod an irreducible poly (lists of ints)
def find_irreducible(p, k):
    x = sympy.symbols('x')
    for coeffs in itertools.product(range(p), repeat=k):
        poly = sympy.Poly([1] + list(coeffs), x, modulus=p)
        if poly.is_irreducible: return [1] + list(coeffs)
def field_ops(p, k, modpoly):
    # elements: tuples of length k (coeffs of 1, x, ..., x^{k-1})
    def add(a, b): return tuple((x + y) % p for x, y in zip(a, b))
    def mul(a, b):
        prod = [0] * (2 * k - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    prod[i + j] = (prod[i + j] + x * y) % p
        for i in range(2 * k - 2, k - 1, -1):
            c = prod[i]
            if c:
                for j in range(k + 1):
                    prod[i - k + j] = (prod[i - k + j] - c * modpoly[::-1][j]) % p  # modpoly monic: x^k + ... ; modpoly[::-1] gives [c0..ck]
        return tuple(prod[:k])
    return add, mul
def count_points(coeffs, p, k):
    """number of F_{p^k}-points of the smooth projective curve w^3 = f(t) (f of degree 3 or 6, monic, squarefree).
    Affine points via cube-root counting: for q = 1 mod 3, #{w : w^3 = c} = 3 if c is a nonzero cube, 1 if c = 0, 0 otherwise;
    for q = 2 mod 3 cubing is a bijection: exactly 1 solution.  Points at infinity: deg f = 6: three (w/t^2)^3 = 1 -> 3 if q = 1 mod 3 else 1;
    deg f = 3: w^3 = v^3 + ... : (w/v)^3 = 1 -> same count."""
    q = p ** k
    if k == 1:
        elems = [(i,) for i in range(p)]; add = lambda a, b: ((a[0] + b[0]) % p,); mul = lambda a, b: ((a[0] * b[0]) % p,)
    else:
        modpoly = find_irreducible(p, k)
        add, mul = field_ops(p, k, modpoly)
        elems = list(itertools.product(range(p), repeat=k))
    zero = tuple([0] * k)
    # cube map image counts
    cubes = {}
    for w in elems:
        c = mul(mul(w, w), w); cubes[c] = cubes.get(c, 0) + 1
    # evaluate f at t by Horner
    def ev(t):
        acc = zero
        for c in reversed(coeffs):
            acc = add(mul(acc, t), tuple([c % p] + [0] * (k - 1)))
        return acc
    n = 0
    for t in elems:
        n += cubes.get(ev(t), 0)
    n += 3 if q % 3 == 1 else 1
    return n
t0 = time.time()
res = {}
for k in (1, 2, 3):
    if p ** k > 60000: break
    nT = count_points(bco, p, k); nE = count_points(gco, p, k)
    res[k] = (nT, nE)
    print(f'p={p} k={k}: #T_B={nT} #E\'={nE}', flush=True)
q = p
# traces: #C(F_{q^k}) = q^k + 1 - s_k(C); Prym traces s_k(P) = s_k(T_B) - s_k(E')
s = {k: (p ** k + 1 - res[k][0]) - (p ** k + 1 - res[k][1]) for k in res}
s1, s2, s3 = s[1], s[2], s[3]
e1 = s1; e2 = (s1 ** 2 - s2) // 2; e3 = (s1 ** 3 - 3 * s1 * s2 + 2 * s3) // 6
assert (s1 ** 2 - s2) % 2 == 0 and (s1 ** 3 - 3 * s1 * s2 + 2 * s3) % 6 == 0
c1, c2, c3 = -e1, e2, -e3
T = sympy.symbols('T')
L = 1 + c1 * T + c2 * T ** 2 + c3 * T ** 3 + p * c2 * T ** 4 + p ** 2 * c1 * T ** 5 + p ** 3 * T ** 6
fac = sympy.factor_list(L)
Lpoly = sympy.Poly(L, T)
# root-of-unity ratio test: resultant of L(x) and L(zx) wrt x against cyclotomic polynomials
irreducible = len(fac[1]) == 1 and fac[1][0][1] == 1
rootsq = [complex(r) for r in sympy.Poly(L, T).nroots(n=30)]
import cmath
absok = all(abs(abs(r) - p ** -0.5) < 1e-9 for r in rootsq)
# check ratios are not roots of unity of order <= 2000 (numerically)
ratios_ok = True
for i in range(6):
    for j in range(6):
        if i == j: continue
        r = rootsq[i] / rootsq[j]
        ang = cmath.phase(r) / (2 * cmath.pi)
        for m in range(1, 2001):
            if abs(ang * m - round(ang * m)) < 1e-7: ratios_ok = False; break
        if not ratios_ok: break
    if not ratios_ok: break
# also E' traces and check j=0 shape: a^2 - 4p = -3 m^2
aE = p + 1 - res[1][1]
disc = aE ** 2 - 4 * p
m2 = -disc // 3 if disc % 3 == 0 else None
ent = dict(exp='a5-prym-independent', p=p, counts=res, prym_traces=s, L_P=str(sympy.expand(L)), factor_over_Q=str(fac),
           irreducible_over_Q=irreducible, roots_abs_ok=absok, no_root_of_unity_ratio_upto_2000=ratios_ok,
           aE=aE, Eprime_j0_shape=(m2 is not None and sympy.sqrt(m2).is_Integer), seconds=round(time.time() - t0, 1))
log(ent)
