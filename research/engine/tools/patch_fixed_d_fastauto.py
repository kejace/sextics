"""Idempotent: solve_fixed_D_fast(method='auto') -- sparse elimination first, compressed rows as fallback."""
import sys
P = '/home/kejace/devel/sextics-work/hodge/fixed_d.py'
src = open(P).read()
if "method='auto'" in src:
    print('fixed_d fastauto: already patched'); sys.exit(0)
old = '''def solve_fixed_D_fast(F, A, q, D, p=None, verbose=False, want_exact=True):
    """Compressed modular consistency test, then exact solve if consistent."""
    p = p or core.PRIMES[0]
    S = FixedDSystem(F, A, q, D, p, verbose=verbose)
    r = solve_compressed(S, verbose=verbose)
    out = {'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
           'rows': S.nrows, 'cols': S.ncols, 'd': S.d}'''
new = '''def solve_fixed_D_fast(F, A, q, D, p=None, verbose=False, want_exact=True, method='auto', sparse_kw=None):
    """Modular consistency test, then exact solve if consistent.

    method: 'auto' (sparse elimination with a fill limit of 25x the initial nnz, the dense finish absorbs the
    rest; compressed rows only if the sparse path raises), 'sparse', or 'compressed' (the pre-v3 behaviour).
    The result carries 'solver' (which path produced the modular answer) and, for sparse, its statistics."""
    p = p or core.PRIMES[0]
    S = FixedDSystem(F, A, q, D, p, verbose=verbose)
    r = None; solver = None
    if method in ('auto', 'sparse'):
        try:
            kw = dict(sparse_kw or {})
            if method == 'auto':
                kw.setdefault('max_nnz', 25 * max(1, sum(len(c) for c in S.colpolys)))
            r = solve_sparse(S, verbose=verbose, **kw); solver = 'sparse'
            if r['consistent'] is None:
                raise RuntimeError(r.get('note', 'sparse solve unverified'))
        except Exception as ex:
            if method == 'sparse':
                raise
            if verbose:
                print(f'[fast] sparse path failed ({ex.__class__.__name__}: {ex}); using compressed rows', file=sys.stderr)
            r = None
    if r is None:
        r = solve_compressed(S, verbose=verbose); solver = 'compressed'
    out = {'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
           'rows': S.nrows, 'cols': S.ncols, 'd': S.d, 'solver': solver}
    if solver == 'sparse':
        out['sparse_stats'] = {k: r.get(k) for k in ('status_sparse', 'npiv_sparse', 'active_rows', 'active_cols',
                                                     'active_nnz', 'fill_events', 't_sparse', 't_dense', 't_total')}'''
assert old in src, 'solve_fixed_D_fast not found'
src = src.replace(old, new)
open(P, 'w').write(src); print('fixed_d fastauto: patched')
