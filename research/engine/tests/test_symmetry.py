"""symmetry.py against accepted witnesses: group sizes, characters, residue character of A,
and the sector-restricted solver must reproduce an exact certificate for the accepted D."""
from common import *
import time
from hodge import core, fixed_d, symmetry

def test_group_fermat_like():
    F, A, q, _ = core.load_instance('fermat-pair-rational-d4-strict')
    G = symmetry.diagonal_group(F)
    assert len(G) == 6 ** 5 and G.order == 6, (len(G), G.order)
    for w in G.elements[:200]:
        assert symmetry.character(F, w, G.order) == 0
    # the residue character of A is 0 on a subgroup (never empty: identity)
    inv = symmetry.residue_invariant_elements(F, A, q, G)
    assert (0,) * 6 in inv and len(G) % len(inv) == 0
    print(f'  fermat-pair: |G|={len(G)} order={G.order} residue-invariant subgroup {len(inv)}')

def test_group_general():
    for iid in ['collision-h3', 'hesse-triangle', 'klein-hesse', 'contact46006']:
        F, A, q, _ = core.load_instance(iid)
        G = symmetry.diagonal_group(F)
        n = G.order
        for w in G.elements:
            assert symmetry.character(F, w, n) is not None
        # closure under addition for a few random pairs
        import random
        S = set(G.elements)
        for _ in range(50):
            a, b = random.choice(G.elements), random.choice(G.elements)
            c = tuple((x + y) % n for x, y in zip(a, b))
            assert c in S
        print(f'  {iid}: |G|={len(G)} order={n} hnf diag={[G.hnf[i][i] for i in range(5)]}')

def test_restricted_solver_reproduces_witness():
    for iid in SMALL + ['hesse-triangle', 'xtype-sparse-rational-1m142']:
        F, A, q, _ = core.load_instance(iid)
        d, D, N, Q = core.load_witness(iid)
        t0 = time.time()
        st, data = symmetry.solve_fixed_D_symmetric(F, A, q, D)
        assert st == 'solved', (iid, st, data)
        ok, info = core.verify_witness(F, A, q, D, data['N'], data['Q'])
        assert ok
        print(f"  {iid}: d={d} restricted cols {data['cols']}/{data['full_cols']} "
              f"(stabilizer {data['subgroup_size']}/{data['group_size']}) solved+verified ({time.time()-t0:.1f}s)")
        # perturbation must be inconsistent (negative control)
        D2 = dict(D); e = next(iter(D2)); D2[e] = D2[e] + 1
        st2, _ = symmetry.solve_fixed_D_symmetric(F, A, q, D2, max_primes=1)
        assert st2 == 'inconsistent', (iid, st2)

def test_sector_matches_full():
    # restricted system: consistent, and its rank is at most the full rank
    F, A, q, _ = core.load_instance('collision-h3')
    d, D, N, Q = core.load_witness('collision-h3')
    S, sub = symmetry.restricted_system(F, A, q, D)
    r = S.solve()
    assert r['consistent']
    full = fixed_d.FixedDSystem(F, A, q, D, core.PRIMES[0]).solve()
    assert full['consistent'] and r['rank'] <= full['rank']
    print(f"  collision-h3: restricted rank {r['rank']} nullity {r['nullity']} vs full rank {full['rank']} nullity {full['nullity']}")

if __name__ == '__main__':
    test_group_fermat_like(); test_group_general(); test_restricted_solver_reproduces_witness(); test_sector_matches_full()
    print('test_symmetry OK')
