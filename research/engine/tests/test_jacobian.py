"""jacobian.py against accepted witnesses: A*D in J(F) for every accepted D (d<=6),
witness D is not in J(F) (a denominator in J would be excluded), F*G is detected as
divisible, random polys are not in J, mult_kernel contains the witness class."""
from common import *
import random, time
from fractions import Fraction
from hodge import core, jacobian

def test_membership_and_divisibility():
    for iid in SMALL + ['hesse-triangle']:
        F, A, q, _ = core.load_instance(iid)
        d, D, N, Q = core.load_witness(iid)
        t0 = time.time()
        ok, info = jacobian.in_jacobian(F, core.mul(A, D))
        assert ok, (iid, info)
        assert not jacobian.in_jacobian(F, D)[0], iid
        G = {m: Fraction(random.randint(-3, 3)) for m in core.monomials(2, 0)}
        ok, Gp = jacobian.divisible_by(F, core.mul(F, G))
        assert ok and Gp == core.reduce_mod(G, core.PRIMES[0]), iid
        assert not jacobian.divisible_by(F, D)[0]
        # J contains m*F_i, and a random degree-6 polynomial is (generically) not in J
        Fi = core.pderiv(F, 3)
        assert jacobian.in_jacobian(F, core.mul(Fi, {core.monomials(1)[0]: Fraction(1)}))[0]
        R = {m: Fraction(random.randint(1, 5)) for m in core.monomials(6, 1)}
        assert not jacobian.in_jacobian(F, R)[0]
        print(f'  {iid}: membership tests ok ({time.time()-t0:.1f}s)')

def test_mult_kernel_contains_witness():
    for iid in SMALL:
        F, A, q, _ = core.load_instance(iid)
        d, D, N, Q = core.load_witness(iid)
        t0 = time.time()
        r = jacobian.mult_kernel(F, A, d, return_reducer=True)
        red = r['reducer']
        p = r['p']
        # witness D modulo J must lie in the span of the kernel basis
        nf = red.normal_form(core.reduce_mod(D, p))
        assert nf, 'witness D reduced to zero mod J'
        from flint import nmod_mat
        sm = red.standard; sidx = {m: i for i, m in enumerate(sm)}
        rows = r['basis'] + [nf]
        flat = [0] * (len(rows) * len(sm))
        for i, v in enumerate(rows):
            for e, c in v.items():
                flat[i * len(sm) + sidx[e]] = c
        M = nmod_mat(len(rows), len(sm), flat, p)
        assert M.rank() == len(r['basis']), (iid, M.rank(), len(r['basis']))
        print(f"  {iid}: d={d} ker[A] dim={r['dim']} (R/J)_d={r['ring_dim']} contains witness ({time.time()-t0:.1f}s)")

def test_normal_form_idempotent():
    F, A, q, _ = core.load_instance('collision-h3')
    P = core.mul(A, {core.monomials(2, 0)[3]: Fraction(1)})  # degree 8 odd
    nf = jacobian.normal_form(F, P)
    nf2 = jacobian.normal_form(F, nf)
    assert nf == nf2
    assert jacobian.in_jacobian(F, core.sub(core.reduce_mod(P, core.PRIMES[0]), nf, core.PRIMES[0]))[0]
    print('  normal form idempotent and P - nf(P) in J')

if __name__ == '__main__':
    test_membership_and_divisibility(); test_mult_kernel_contains_witness(); test_normal_form_idempotent()
    print('test_jacobian OK')
