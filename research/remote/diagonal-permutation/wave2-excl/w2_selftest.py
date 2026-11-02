"""Self-test of w2lib: positive control (accepted witness), negative control (wave-1 inconsistent D),
relaxation lemma check, and d=5 relaxation.  Usage: W2IID=<iid> python3 w2_selftest.py <mode>
modes: control (needs W2IID set to an accepted instance with witness.json), lane (diagonal-permutation)."""
import sys, os, time, json
mode = sys.argv[1] if len(sys.argv) > 1 else 'lane'
import w2lib as W
from hodge import core, fixed_d
from hodge.core import monomials, mul, add, sub, scale, pderiv, reduce_mod, degree
import numpy as np

p = W.P1
if mode == 'control':
    d, D, N, Q = core.load_witness(W.IID)
    print('instance', W.IID, 'witness d =', d)
    # (1) E1 via solve_E1 and the lemma: Psi(D;u) in J (degree 2d)
    u = W.solve_E1(D, p)
    assert u is not None, 'E1 failed on an accepted witness?!'
    # check u(F) = A D mod p
    Fp = reduce_mod(W.F, p); Ap = reduce_mod(W.A, p); Dp = reduce_mod(D, p)
    lhs = W.field_apply(u, Fp, p); rhs = mul(Ap, Dp, p)
    assert sub(lhs, rhs, p) == {}, 'u(F) != A D'
    Ps = W.psi(Dp, u, p)
    par = core.parity_of(Ps)
    red = W.Red(2 * d, p, par, None, verbose=True)
    print('Psi(D;u) in J ?', red.in_J(Ps), '(lemma predicts True for an accepted witness)')
    # also with the witness N as u' = N - (H/6)E: Psi' = D div N - N(D) - D H must be divisible by F
    # (2) reduced E2 system consistent?
    t = time.time()
    S = W.ReducedE2(D, p, u=u, verbose=True)
    r = S.solve(compressed=True, verbose=True)
    print('ReducedE2 consistent (expected True):', r['consistent'], round(time.time() - t, 1), 's')
    # (3) engine full system
    out = fixed_d.solve_fixed_D_fast(W.F, W.A, W.q, D, p=p, want_exact=False)
    print('engine full system consistent (expected True):', out['consistent'], 'rank', out['rank'], 'nullity', out['nullity'])
    # (4) a perturbed D that should fail: D + X0^d
    Dbad = add(D, {tuple([d, 0, 0, 0, 0, 0]): 1})
    ub = W.solve_E1(Dbad, p)
    print('perturbed D passes E1?', ub is not None)
    if ub is not None:
        Sb = W.ReducedE2(Dbad, p, u=ub); rb = Sb.solve(); outb = fixed_d.solve_fixed_D_fast(W.F, W.A, W.q, Dbad, p=p, want_exact=False)
        print('perturbed: reduced consistent', rb['consistent'], 'engine consistent', outb['consistent'])
else:
    # lane: negative control D = F_X0 (d = 5): wave 1 says inconsistent, served only zero class
    Fp = reduce_mod(W.F, p)
    D = pderiv(Fp, 0, p)
    u = W.solve_E1(D, p); assert u is not None
    # u should be A * e_0 up to Koszul: check Psi in J? For D = F_0: Psi = -A F_00 mod J
    Ps = W.psi(D, u, p)
    red = W.Red(10, p, core.parity_of(Ps), None, verbose=True)
    print('Psi(F_X0) in J ?', red.in_J(Ps), '(expected False: no degree-5 certificate)')
    t = time.time()
    S = W.ReducedE2(D, p, u=u, verbose=True); r = S.solve(verbose=True)
    print('ReducedE2(F_X0) consistent:', r['consistent'], round(time.time() - t, 1), 's (wave-1: inconsistent)')
    # (01)-sector version: F_X0 is not a (01)-eigenvector; use F_X0 + F_X1 (invariant) and F_X0 - F_X1 (anti)
    for name, Dd, sg in [('F_X0+F_X1', add(pderiv(Fp, 0, p), pderiv(Fp, 1, p), p), 1), ('F_X0-F_X1', sub(pderiv(Fp, 0, p), pderiv(Fp, 1, p), p), -1)]:
        uu = W.solve_E1(Dd, p)
        S1 = W.ReducedE2(Dd, p, u=uu, verbose=True); r1 = S1.solve(verbose=True)
        S2 = W.ReducedE2(Dd, p, u=uu, sigma=True, verbose=True); r2 = S2.solve(verbose=True)
        print(name, 'full-sector consistent', r1['consistent'], 'rank', r1['rank'], '| (01)-sector consistent', r2['consistent'], 'rank', r2['rank'], 'cols', r2['cols'])
    # relaxation at d = 5: 6 unknowns, must have only the trivial solution
    t = time.time()
    Rx = W.Relaxation(5, [0, 1], p, verbose=True)
    E = Rx.independent_equations()
    print('d=5 relaxation: unknowns', Rx.n, 'independent quadratic equations', E.shape[0], round(time.time() - t, 1), 's')
    Rx.write_msolve('/tmp/hodge/results/diagonal-permutation/wave2-excl/relax_d5_mixed.ms')
    # sanity: the Euler field is a solution of the relaxation at d=6 (even sector, no complement)
    Rx6 = W.Relaxation(6, [0], p, euler_complement=False, verbose=True)
    x = np.zeros(Rx6.n, dtype=np.int64)
    for k, (kind, data) in enumerate(Rx6.unk):
        i, m = data
        if m[i] == 1 and sum(m) == 1:
            x[k] = 1
    res = Rx6.evaluate(x)
    print('d=6 even relaxation: unknowns', Rx6.n, 'raw equations', Rx6.neq_raw, '; Euler field residual zero?', not np.any(res))
    y = np.zeros(Rx6.n, dtype=np.int64); y[0] = 1
    print('random unit vector residual nonzero?', bool(np.any(Rx6.evaluate(y))))
