"""E1-level served space T1(D) = {A' in S_6^- : A' D in J_12} for pure-parity degree-6 denominators D (mod p).
served_Q(D) (targets with a certificate) is contained in T1_Q(D), and dim T1_Q <= dim T1_p (rank inequality, with the
Hilbert check for J_12).  J_6^- (dim 18) is always served (N = a D with a linear traceless; Q = 0), so T1 = J_6^- means
D serves no primitive target at all; by the book's q=1 sum principle served(D1 D2) = served(D1)+served(D2) (and
served(D^k) = served(D)), products of such D's serve nothing primitive either.
Usage: python3 w2x_t1.py <family-set>   (same sets as w2x_e1screen; only pure-parity candidates of degree 6)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave2-excl')
from w2x_lib import *
from w2x_quad import rref_J, normal_forms
from w2x_e1screen import families
from hodge.core import monomials, mul, degree, parity_of

p = P
if __name__ == '__main__':
    which = sys.argv[1]
    F, A, q, rec = load(); Fp = core.reduce_mod(F, p)
    cands = [(nm, fam, D) for nm, fam, D in families(which) if degree(D) == 6 and parity_of(D) is not None]
    print(f'{len(cands)} pure-parity degree-6 candidates', flush=True)
    odd6 = monomials(6, 1)
    red = {}
    t0 = time.time()
    for par in (0, 1):
        red[par] = rref_J(Fp, 12, par)
        assert red[par][4] == expected_rank_J(12, par), 'Hilbert check failed'
    print(f'J_12 rref both parities ({time.time()-t0:.0f}s)', flush=True)
    # J_6^- basis to confirm containment
    out = []
    for nm, fam, D in cands:
        t1 = time.time()
        Dp = core.reduce_mod(D, p); parD = parity_of(Dp); par = (1 + parD) % 2
        mons, idx, Rr, pivots, rk = red[par]
        polys = [mul({m: 1}, Dp, p) for m in odd6]
        blocks = normal_forms(polys, mons, idx, Rr, pivots)
        # stack blocks -> matrix (226 x dim S_12^par); kernel of its transpose action = T1
        U = blocks[0]
        for b in blocks[1:]: U = U.concat(b) if False else U  # single block (226 rows < chunk)
        rank = U.rank()
        dimT1 = len(odd6) - rank
        rec_ = {'exp': 't1', 'name': nm, 'family': fam, 'd': 6, 'parD': parD, 'p': p, 'dim T1': dimT1, 'dim J_6^-': 18,
                'T1_trivial': dimT1 == 18, 'seconds': round(time.time() - t1, 1)}
        out.append(rec_); log(rec_)
        print(json.dumps(rec_), flush=True)
    with open(f'{PK}/t1_{which}.json', 'w') as f: json.dump(out, f, indent=1)
    nontriv = [o['name'] for o in out if not o['T1_trivial']]
    summary = {'exp': 't1-summary', 'set': which, 'n': len(out), 'nontrivial_T1': nontriv, 'p': p}
    log(summary); print(json.dumps(summary))
