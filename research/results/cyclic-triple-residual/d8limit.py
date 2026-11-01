"""D8(t) for the family F_t = Fermat + t*couplings: exact kernel element at several t, canonical representative
with respect to the Fermat-ring monomial basis (exponents <= 4 in the sector), to read off the t -> 0 limit."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from sect import *
WD = '/tmp/hodge/results/cyclic-triple-residual'
A = core.parse('X0**3*Y0**3')
coup = core.parse('X0**3*X1*Y0*Y1 + X0**3*X2**3 + X0**2*Y0*Y1*Y2**2 + X0*X1*Y1**3*Y2 + X0*X2**3*Y1*Y2 + X1*X2**3*Y0*Y1')
Ferm = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
primes = PRIMES[:6]
def D8_of(tv):
    Ft = core.add(Ferm, core.scale(coup, tv))
    Fz = integral(Ft); G = semi_invariant_subgroup(elems(Fz), [integral(A)]); gens = generators(G); chi = make_chi(gens)
    c = (0, 4)
    Dm = [m for m in mons(8) if chi(m) == c]; idx = {m: i for i, m in enumerate(Dm)}
    # order columns so that non-Fermat-basis monomials (some exponent >= 5) come first -> pivots of J land there, and the
    # canonical representative is expressed in Fermat-basis monomials (all exponents <= 4) as far as possible
    order = sorted(range(len(Dm)), key=lambda i: (0 if max(Dm[i]) >= 5 else 1, Dm[i]))
    Dm_o = [Dm[i] for i in order]; idx_o = {m: i for i, m in enumerate(Dm_o)}
    dF = [ideriv(Fz, i) for i in range(6)]
    vecs_by_p = {}
    for p in primes:
        res, _ = e1_kernel(Ft, A, 8, p=p, G=G)
        r = res[c]; assert r['ker_mod_J'] == 1, r['ker_mod_J']
        Jrows = []
        for i in range(6):
            ci = cadd(c, chi(E6[i]))
            for m in mons(3):
                if chi(m) != ci: continue
                row = [0]*len(Dm_o)
                for e, v in dF[i].items():
                    k = tuple(a+b for a, b in zip(m, e)); row[idx_o[k]] = (row[idx_o[k]] + v) % p
                Jrows.append(row)
        J = nmod_mat(len(Jrows), len(Dm_o), [x for r_ in Jrows for x in r_], p); Jr, rk = J.rref()
        piv = []
        for i in range(rk):
            j = 0
            while int(Jr[i, j]) == 0: j += 1
            piv.append(j)
        best = None
        for vec in r['vectors']:
            v = [0]*len(Dm_o)
            for m, x in vec.items(): v[idx_o[m]] = x % p
            for i, j in enumerate(piv):
                if v[j]:
                    f = v[j] * pow(int(Jr[i, j]), -1, p) % p
                    for k in range(len(Dm_o)):
                        if int(Jr[i, k]): v[k] = (v[k] - f*int(Jr[i, k])) % p
            if any(v): best = v; break
        # normalise: coefficient of the LAST nonzero monomial (in Fermat-basis order) = 1
        j0 = max(k for k in range(len(Dm_o)) if best[k]); inv = pow(best[j0], -1, p)
        vecs_by_p[p] = [(x*inv) % p for x in best]
    acc, m = vecs_by_p[primes[0]], primes[0]
    for p in primes[1:]:
        acc = [core.crt_pair(a, m, b, p)[0] for a, b in zip(acc, vecs_by_p[p])]; m *= p
    rat = [core.ratrecon(x, m) for x in acc]
    if any(x is None for x in rat): return None, Dm_o
    return {Dm_o[i]: rat[i] for i in range(len(Dm_o)) if rat[i]}, Dm_o
for tv in [Fraction(5, 7), Fraction(1, 10), Fraction(1, 100)]:
    t0 = time.time()
    D, Dm_o = D8_of(tv)
    if D is None:
        print(f't={tv}: rational reconstruction failed', flush=True); continue
    nonF = {m: c for m, c in D.items() if max(m) >= 5}
    print(f't={tv} ({time.time()-t0:.0f}s): {len(D)} terms, {len(nonF)} with an exponent >= 5; normalised representative:', flush=True)
    for m, c in sorted(D.items(), key=lambda kv: -abs(float(kv[1]))):
        print(f'    {str(c):>22s}  {core.to_str({m: 1})}', flush=True)
    json.dump({'t': str(tv), 'D': core.to_str(D)}, open(f'{WD}/d8_t_{tv.numerator}_{tv.denominator}.json', 'w'))
