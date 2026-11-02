"""Exact rational representatives of the invariant E1-kernel classes via CRT over several primes."""
import sys, time, json, pickle, itertools
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power

def exact_kernel(d, pD, nprimes=5):
    primes = core.PRIMES[:nprimes]
    reps_by_prime = []
    for p in primes:
        info, polys = e1_kernel_inv(F, A, q, d, pD, p=p)
        reps_by_prime.append([core.reduce_mod(P, p) for P in polys])
    k = len(reps_by_prime[0])
    assert all(len(r) == k for r in reps_by_prime), 'kernel dimension differs across primes'
    exact = []
    for i in range(k):
        keys = set()
        for r in reps_by_prime: keys |= set(r[i].keys())
        P = {}
        for e in keys:
            acc, m = None, None
            for r, p in zip(reps_by_prime, primes):
                v = r[i].get(e, 0)
                if acc is None: acc, m = v, p
                else: acc, m = core.crt_pair(acc, m, v, p)
            fr = core.ratrecon(acc, m)
            if fr is None:
                return None, 'ratrecon failed'
            if fr: P[e] = fr
        exact.append(P)
    return exact, 'ok'

if __name__ == '__main__':
    degs = [(8, 0), (9, 0), (9, 1)] if len(sys.argv) < 2 else [tuple(map(int, a.split(','))) for a in sys.argv[1:]]
    e = {('X', k): e_k(X, k) for k in (1, 2, 3)}; e.update({('Y', k): e_k(Y, k) for k in (1, 2, 3)})
    for d, pD in degs:
        t0 = time.time()
        exact, msg = exact_kernel(d, pD)
        print(f'=== d={d} pD={pD}: {msg} ({time.time()-t0:.1f}s)', flush=True)
        if exact is None: continue
        pickle.dump(exact, open(f'{RES}/exactK_d{d}_p{pD}.pkl', 'wb'))
        # e-monomial coordinates (exact) via fmpq solve
        from flint import fmpq_mat, fmpq
        emons = []
        for a in itertools.product(range(d + 1), repeat=6):
            w = a[0] + 2*a[1] + 3*a[2] + a[3] + 2*a[4] + 3*a[5]
            if w != d or (a[3] + 2*a[4] + 3*a[5]) % 2 != pD: continue
            P = {(0,)*6: Fraction(1)}
            for kk, (blk, j) in enumerate([('X',1),('X',2),('X',3),('Y',1),('Y',2),('Y',3)]):
                if a[kk]: P = mul(P, pw(e[(blk, j)], a[kk]))
            emons.append((a, P))
        Dreps = orbit_reps(d, pD, canon_codes)
        n = len(Dreps); m = len(emons)
        M = fmpq_mat(n, m)
        for c, (_, P) in enumerate(emons):
            for r, rep in enumerate(Dreps):
                v = P.get(rep, 0)
                if v: M[r, c] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
        for i, P in enumerate(exact):
            ok, info = fixed_d.e1_solvable(F, A, q, P)
            b = fmpq_mat(n, 1)
            for r, rep in enumerate(Dreps):
                v = P.get(rep, 0)
                if v: b[r, 0] = fmpq(Fraction(v).numerator, Fraction(v).denominator)
            try:
                x = M.solve(b)
                terms = [(emons[j][0], Fraction(int(x[j,0].p), int(x[j,0].q))) for j in range(m) if x[j,0] != 0]
            except Exception as ex:
                terms = [('solve failed', str(ex))]
            print(f'--- rep {i}: {len(P)} monomials, E1 exact(mod p) = {ok}; max |num| = {max(abs(Fraction(v).numerator) for v in P.values())}, lcm den = {max(Fraction(v).denominator for v in P.values())}')
            for a, v in terms:
                if isinstance(a, tuple):
                    print('    ', v, ' * e1X^%d e2X^%d e3X^%d e1Y^%d e2Y^%d e3Y^%d' % a)
                else: print('    ', a, v)
