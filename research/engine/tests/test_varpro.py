"""varpro.py: gradient check, recovery of accepted denominators from monomial families,
spurious (F | D) detection, exact confirmation."""
from common import *
import time
import numpy as np
from fractions import Fraction
from hodge import core, varpro, fixed_d

def test_gradient():
    F, A, q, _ = core.load_instance('collision-h3')
    fam = varpro.VarProFamily(F, A, q, varpro.monomial_family(2, 0))
    rng = np.random.default_rng(1)
    u = rng.standard_normal(fam.dimP); u /= np.linalg.norm(u)
    f0, g = fam.value_and_grad(u)
    h = 1e-5; gf = np.zeros_like(u)
    for k in range(len(u)):
        e = np.zeros_like(u); e[k] = h
        gf[k] = (fam.value_and_grad(u + e)[0] - fam.value_and_grad(u - e)[0]) / (2 * h)
    err = np.max(np.abs(g - gf)) / np.max(np.abs(g))
    assert err < 1e-6, err
    print(f'  gradient check rel err {err:.2e}')

def test_recover_d2_and_d4():
    for iid, extra in [('collision-h3', 0), ('fermat-pair-rational-d4-strict', 6), ('three-block-rational-d4-strict', 6)]:
        F, A, q, _ = core.load_instance(iid)
        d, D, N, Q = core.load_witness(iid)
        if d == 2:
            basis = varpro.monomial_family(2, 0)
        else:
            # witness monomials plus `extra` unrelated even monomials
            ms = list(D)
            others = [m for m in core.monomials(d, 0) if m not in D][:extra]
            basis = [{m: Fraction(1)} for m in ms + others]
        t0 = time.time()
        res = varpro.search_and_confirm(F, A, q, basis, n_restarts=12, seed=0)
        assert res['best_rel'] < 1e-10, (iid, res['best_rel'])
        assert res['confirmed'], (iid, 'no exact confirmation', res['tried'])
        c = res['confirmed'][0]
        ok, _ = core.verify_witness(F, A, q, c['D'], c['N'], c['Q'])
        assert ok
        print(f"  {iid}: K={len(basis)} best rel {res['best_rel']:.1e} confirmed D = {core.to_str(c['D'])[:70]} ({time.time()-t0:.1f}s)")

def test_spurious_detection():
    F, A, q, _ = core.load_instance('hesse-triangle')
    G = {core.monomials(0)[0]: Fraction(1)}
    basis = [core.mul(F, G), {core.monomials(6, 0)[5]: Fraction(1)}, {core.monomials(6, 0)[9]: Fraction(1)}]
    fam = varpro.VarProFamily(F, A, q, basis, exclude_F_multiples=True)
    assert fam.excluded_dim == 1 and fam.dimP == 2
    spur, rel = fam.is_spurious(np.array([1.0, 0.0, 0.0]))
    assert spur and rel < 1e-8
    spur2, rel2 = fam.is_spurious(np.array([0.0, 1.0, 1.0]))
    assert not spur2
    # the F-multiple direction is orthogonal to the reduced parameter space
    assert np.allclose(fam.P.T @ np.array([1.0, 0, 0]), 0)
    print(f'  spurious detection ok (excluded dim {fam.excluded_dim}, div rel {rel:.1e} vs {rel2:.1e})')

def test_reconstruct():
    F, A, q, _ = core.load_instance('collision-h3')
    fam = varpro.VarProFamily(F, A, q, varpro.monomial_family(2, 0))
    a = np.array([1 / 3, -2 / 7, 0.5, 1.0, 0, 0, 0, 0, 0, 0, 0, 0])
    cands = fam.reconstruct(a, max_den=50)
    assert any(tuple(c[:4]) == (Fraction(1, 3), Fraction(-2, 7), Fraction(1, 2), Fraction(1)) for c in cands)
    print('  rational reconstruction ok')

if __name__ == '__main__':
    test_gradient(); test_spurious_detection(); test_reconstruct(); test_recover_d2_and_d4()
    print('test_varpro OK')
