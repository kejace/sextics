"""Point counts on C_B: w^3 = b(t), b(t) = 1 + 15/7 t^2 + 10/7 t^4 + t^6 (genus 4), its elliptic quotient
E: w^3 = 1 + 15/7 v + 10/7 v^2 + v^3 (v = t^2), and the Prym P = ker(J(C_B) -> E) (dim 3).
We compute L-polynomials over F_p (p = 1 mod 3, p != 7) from #C(F_{p^k}), k=1..4, and factor L_P = L_C / L_E over Q."""
import sys
from sympy import symbols, Poly, factor_list, ZZ, GF
from sympy import factorint
T = symbols('T')
def gf_count_curve(p, k, which):
    """#points of the smooth projective curve over F_{p^k}; 'C': w^3=b(t); 'E': w^3=e(v)."""
    # Build F_{p^k} via a Conway-free approach: use sympy's GF? too slow; implement F_{p^k} as polynomials mod irreducible.
    import itertools
    # find an irreducible monic poly of degree k over F_p
    from sympy import Poly as SP
    x = symbols('x')
    if k == 1:
        elems = list(range(p))
        def mul(a, b): return a * b % p
        def add(a, b): return (a + b) % p
        def pw(a, n): return pow(a, n, p)
        zero, one = 0, 1
        inv7 = pow(7, -1, p)
        def const(c): return c % p
    else:
        # represent elements as tuples of length k; find irreducible
        import random
        random.seed(1)
        while True:
            coeffs = [random.randrange(p) for _ in range(k)] + [1]
            f = SP(coeffs[::-1], x, modulus=p)
            if f.is_irreducible: break
        mod = [c % p for c in coeffs]  # low->high, monic degree k
        def add(a, b): return tuple((u + v) % p for u, v in zip(a, b))
        def mul(a, b):
            res = [0] * (2 * k - 1)
            for i, u in enumerate(a):
                if u == 0: continue
                for j, v in enumerate(b):
                    res[i + j] = (res[i + j] + u * v) % p
            for d in range(2 * k - 2, k - 1, -1):
                c = res[d]
                if c:
                    for i in range(k):
                        res[d - k + i] = (res[d - k + i] - c * mod[i]) % p
                    res[d] = 0
            return tuple(res[:k])
        def pw(a, n):
            r = tuple([1] + [0] * (k - 1)); b = a
            while n:
                if n & 1: r = mul(r, b)
                b = mul(b, b); n >>= 1
            return r
        zero = tuple([0] * k); one = tuple([1] + [0] * (k - 1))
        def const(c): return tuple([c % p] + [0] * (k - 1))
        elems = list(itertools.product(range(p), repeat=k))
        inv7 = pow(7, -1, p)
    q = p ** k
    # number of cube roots: for y in F_q, #{w: w^3 = y} = 1 if y=0, else 3 if y is a cube (q=1 mod 3) else 0
    # count cubes via the map; precompute set of cubes
    cubes = {}
    for w in elems:
        c = pw(w, 3)
        cubes[c] = cubes.get(c, 0) + 1
    c15 = const(15 * inv7); c10 = const(10 * inv7)
    total = 0
    for t in elems:
        if which == 'C':
            t2 = mul(t, t); t4 = mul(t2, t2); t6 = mul(t4, t2)
            y = add(add(add(one, mul(c15, t2)), mul(c10, t4)), t6)
        else:
            v = t; v2 = mul(v, v); v3 = mul(v2, v)
            y = add(add(add(one, mul(c15, v)), mul(c10, v2)), v3)
        total += cubes.get(y, 0)
    # points at infinity: C: w^3 = b(t), deg 6, in P(1,1,2): at t=inf, w^3 = 1 (leading coeff 1): 3 points (if q=1 mod 3).
    # E: w^3 = v^3 + ..., at infinity: one point (smooth cubic model: w^3 = v^3 + ... has points at infinity where w^3 = v^3: 3 points w = zeta v) -> 3 points.
    total += 3
    return total
def lpoly(counts, p, g):
    # from N_k = q^k + 1 - sum alpha_i^k, k=1..g -> power sums -> coefficients of L(T) = prod(1 - alpha_i T)
    s = [p ** k + 1 - counts[k - 1] for k in range(1, g + 1)]  # power sums of alphas
    e = [1]
    for k in range(1, g + 1):
        ek = sum((-1) ** (i - 1) * e[k - i] * s[i - 1] for i in range(1, k + 1)) // k
        e.append(ek)
    # L(T) = sum_{i} (-1)^i e_i T^i for i<=g, and L_{2g-i} = p^{g-i} L_i
    coeffs = [(-1) ** i * e[i] for i in range(g + 1)]
    full = coeffs + [p ** (i - g) * coeffs[2 * g - i] for i in range(g + 1, 2 * g + 1)]
    return Poly(sum(c * T ** i for i, c in enumerate(full)), T)
for p in [13, 19, 31, 37, 43]:
    if p % 3 != 1 or p == 7: continue
    NC = [gf_count_curve(p, k, 'C') for k in range(1, 5)]
    NE = [gf_count_curve(p, k, 'E') for k in range(1, 2)]
    LC = lpoly(NC, p, 4); LE = lpoly(NE, p, 1)
    LP, rem = LC.div(LE)
    print(f'p={p}: #C={NC} #E={NE}')
    print(f'   L_E = {LE.as_expr()}')
    print(f'   L_P = {LP.as_expr()}   (rem {rem.as_expr()})')
    print(f'   factor L_P over Q: {factor_list(LP.as_expr())}')
    sys.stdout.flush()
