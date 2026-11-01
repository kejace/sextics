"""L-polynomials of C_B: w^3 = b(t) (genus 4), E: w^3 = b2(v) (v = t^2), Prym P = J(C_B)/E (dim 3), over F_p, p = 1 mod 3.
Uses #C(F_{p^k}), k = 1..3 and L_C = L_E * L_P with L_P of degree 6 (3 free coefficients)."""
import sys, itertools, random
from sympy import symbols, Poly, factor_list
T = symbols('T')
def field(p, k):
    if k == 1:
        return list(range(p)), (lambda a, b: a * b % p), (lambda a, b: (a + b) % p), (lambda c: c % p)
    random.seed(1)
    from sympy import Poly as SP
    x = symbols('x')
    while True:
        coeffs = [random.randrange(p) for _ in range(k)] + [1]
        if SP(coeffs[::-1], x, modulus=p).is_irreducible: break
    mod = [c % p for c in coeffs]
    def add(a, b): return tuple((u + v) % p for u, v in zip(a, b))
    def mul(a, b):
        res = [0] * (2 * k - 1)
        for i, u in enumerate(a):
            if u:
                for j, v in enumerate(b):
                    res[i + j] = (res[i + j] + u * v) % p
        for d in range(2 * k - 2, k - 1, -1):
            c = res[d]
            if c:
                for i in range(k): res[d - k + i] = (res[d - k + i] - c * mod[i]) % p
                res[d] = 0
        return tuple(res[:k])
    def const(c): return tuple([c % p] + [0] * (k - 1))
    return list(itertools.product(range(p), repeat=k)), mul, add, const
def count(p, k, which):
    elems, mul, add, const = field(p, k)
    cubes = {}
    for w in elems:
        c = mul(mul(w, w), w); cubes[c] = cubes.get(c, 0) + 1
    i7 = pow(7, -1, p); c15 = const(15 * i7); c10 = const(10 * i7); one = const(1)
    tot = 0
    for t in elems:
        if which == 'C':
            t2 = mul(t, t); t4 = mul(t2, t2); t6 = mul(t4, t2)
            y = add(add(add(one, mul(c15, t2)), mul(c10, t4)), t6)
        else:
            v2 = mul(t, t); v3 = mul(v2, t)
            y = add(add(add(one, mul(c15, t)), mul(c10, v2)), v3)
        tot += cubes.get(y, 0)
    return tot + 3   # three points at infinity (leading coefficient 1 is a cube; q = 1 mod 3)
def lpoly_from_powersums(s, p, g):
    e = [1]
    for k in range(1, g + 1):
        e.append(sum((-1) ** (i - 1) * e[k - i] * s[i - 1] for i in range(1, k + 1)) // k)
    c = [(-1) ** i * e[i] for i in range(g + 1)]
    full = c + [p ** (i - g) * c[2 * g - i] for i in range(g + 1, 2 * g + 1)]
    return Poly(sum(ci * T ** i for i, ci in enumerate(full)), T)
for p in [13, 19, 31, 37, 43, 61, 67]:
    NC = [count(p, k, 'C') for k in (1, 2, 3)]
    NE = count(p, 1, 'E')
    aE = p + 1 - NE
    LE = Poly(1 - aE * T + p * T ** 2, T)
    # power sums of Frobenius eigenvalues on C: s_k = p^k + 1 - N_k ; on E: aE-powers via Newton; on P: difference
    sC = [p ** k + 1 - NC[k - 1] for k in (1, 2, 3)]
    # E eigenvalues power sums
    sE = [aE, aE ** 2 - 2 * p, aE ** 3 - 3 * p * aE]
    sP = [sC[i] - sE[i] for i in range(3)]
    LP = lpoly_from_powersums(sP, p, 3)
    LC = lpoly_from_powersums(sC + [0], p, 4) if False else None
    fl = factor_list(LP.as_expr())
    print(f'p={p}: #C(F_p,F_p^2,F_p^3)={NC}, #E(F_p)={NE}, a_E={aE}')
    print(f'   L_P = {LP.as_expr()}')
    print(f'   factorization over Q: {fl}', flush=True)
