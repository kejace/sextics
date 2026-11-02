"""Rigorous fixed-D exclusion (q = 1) via the trivial-kernel rank criterion.

For a pure-parity D of degree d the fixed-D system  M x = b  (x = coefficients of N_0..N_5, Q) has the explicit
rational homogeneous solutions  N = F V + G z,  Q = D div V - V(D)   (V in S_{d-4}^6 with the numerator parities,
G in S_d with the parity of Q), whose span has dimension  triv(d) = 6 s_{d-4} + s_d^{par Q} - s_{d-6}^{par Q}
(kernel of (V,G) -> N is {(M z, -F M)}; F irreducible).  Hence over Q:  nullity_Q >= triv.  Mod p (p not dividing
the data denominators):  rank_p <= rank_Q, so nullity_p >= nullity_Q >= triv.  The compressed system (random
row combinations) has rank_c <= rank_p.  THEREFORE: if the compressed nullity C - rank_c equals triv, then
rank_c = rank_p = rank_Q, and an inconsistency of the compressed system (which implies inconsistency of the full
system mod p) implies inconsistency over Q, hence over C (linear system with rational data).  That is a rigorous
exclusion of D (and, at q = 1, of every denominator with the same support, by the prescribed-support theorem).

usage: python3 fixedD_rigor.py <name> <D-source> [--primes k] [--family tag]
   D-source: path to JSON with key 'D' (expression string) or a sympy-style expression."""
import sys, time, json, argparse
from common import *
from hodge import fixed_d
from hodge.core import monomials
ap = argparse.ArgumentParser()
ap.add_argument('name'); ap.add_argument('src'); ap.add_argument('--primes', type=int, default=2)
ap.add_argument('--family', default='rigor'); ap.add_argument('--expr', action='store_true')
ap.add_argument('--k-extra', type=int, default=40)
args = ap.parse_args()
F, A, q, rec = core.load_instance(IID)
if args.expr: D = core.parse(args.src)
else:
    obj = json.load(open(args.src))
    D = core.parse(obj['D']) if isinstance(obj, dict) else core.parse(obj)
d = core.degree(D); pD = core.parity_of(D)
assert pD is not None and core.is_homogeneous(D)
pars = fixed_d.numerator_parities(D, A, F)      # parity of N_i
pQ = core.parity_of(A)
def s(k, par=None): return len(monomials(k, par)) if k >= 0 else 0
triv = sum(s(d - 4, (pars[i] - (0 if i < 3 else 1)) % 2) for i in range(6)) + s(d, pQ) - s(d - 6, pQ)
# V_i parity: N_i = F V_i + G z_i ; parity(V_i) = parity(N_i) ; G parity = parity(N_i) - parity(z_i) = pQ
res = {'exp': 'fixedD-rigor', 'name': args.name, 'family': args.family, 'd': d, 'parityD': pD, 'triv': triv, 'runs': []}
rigorous = None
for p in core.PRIMES[:args.primes]:
    t0 = time.time()
    S = fixed_d.FixedDSystem(F, A, q, D, p, verbose=False)
    r = fixed_d.solve_compressed(S, k_extra=args.k_extra, verbose=False)
    run = {'p': p, 'consistent': r['consistent'], 'rank_compressed': r['rank'], 'nullity_compressed': r['nullity'],
           'rows': S.nrows, 'cols': S.ncols, 'k': r['k'], 'seconds': round(time.time() - t0, 1)}
    run['rank_equality_certified'] = (r['nullity'] == triv)
    run['rigorous_exclusion'] = (not r['consistent']) and run['rank_equality_certified']
    res['runs'].append(run)
    print(run, flush=True)
    if run['rigorous_exclusion']: rigorous = True
    if r['consistent']:
        rigorous = False
        res['note'] = 'CONSISTENT mod p -> attempt exact lift'
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, verbose=True, want_exact=True)
        if out.get('exact'):
            res['note'] += ' SOLVED-EXACT WITNESS'
            json.dump({'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])},
                      open(os.path.join(LANE, f'witness-candidate-{args.name}.json'), 'w'), indent=1)
        break
res['rigorous_exclusion_over_C'] = bool(rigorous)
res['D_head'] = core.to_str(D)[:200]
log(res)
receipt(f'fixedD_{args.name}.json', {**res, 'D': core.to_str(D)})
