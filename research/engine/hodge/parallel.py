"""Multiprocessing over primes: CRT + rational reconstruction of a fixed-D certificate with several primes
solved concurrently (fork), using the sparse solver by default.

    from hodge import parallel
    status, data = parallel.solve_fixed_D_parallel(F, A, q, D, workers=4, method='sparse')

Returns exactly like fixed_d.solve_fixed_D: (status, data) with status in {'inconsistent', 'solved',
'lift-failed'}; data has rank/nullity/rows/cols (+ N, Q when solved).  Each prime runs in its own forked
process (one core each; set FLINT_NUM_THREADS=1 there), batches of `workers` primes are solved, reconstruction
is attempted after every batch.  Pivot structures that differ between primes restart the accumulation (as in
fixed_d.solve_fixed_D).
"""
from __future__ import annotations
import os, sys, time
import multiprocessing as mp
from fractions import Fraction
from . import core, fixed_d
from .core import NV

_CTX = None


def _solve_one(args):
    F, A, q, D, p, method, kw = args
    os.environ['HODGE_QUIET'] = '1'
    try:
        import flint; flint.ctx.threads = 1
    except Exception:
        pass
    S = fixed_d.FixedDSystem(F, A, q, D, p)
    if method == 'dense':
        r = S.solve()
    elif method == 'compressed':
        r = fixed_d.solve_compressed(S)
    else:
        from . import sge
        r = sge.solve_fixed_system(S, **kw)
    return {'p': p, 'consistent': r['consistent'], 'rank': r['rank'], 'nullity': r['nullity'],
            'pivots': r['pivots'], 'x': r['x'], 'rows': S.nrows, 'cols': S.ncols, 'cols_list': S.cols}


def solve_fixed_D_parallel(F, A, q, D, primes=None, max_primes=8, workers=4, method='sparse', verbose=False,
                           sparse_kw=None):
    primes = list(primes or core.PRIMES)[:max_primes]
    d = core.degree(D)
    acc = None; pivots0 = None; last = None
    ctx = mp.get_context('fork')
    done = 0
    with ctx.Pool(min(workers, len(primes))) as pool:
        for start in range(0, len(primes), workers):
            batch = primes[start:start + workers]
            t0 = time.time()
            results = pool.map(_solve_one, [(F, A, q, D, p, method, sparse_kw or {}) for p in batch])
            if verbose:
                print(f'[parallel] batch {batch} solved in {time.time()-t0:.1f}s', file=sys.stderr)
            for r in results:
                done += 1
                base = {'prime': r['p'], 'rank': r['rank'], 'nullity': r['nullity'], 'rows': r['rows'], 'cols': r['cols'], 'd': d}
                if not r['consistent']:
                    pool.terminate()
                    return 'inconsistent', base
                last = base
                if pivots0 is None:
                    pivots0 = r['pivots']
                elif pivots0 != r['pivots']:
                    if verbose:
                        print('[parallel] pivot structure differs between primes; restarting accumulation', file=sys.stderr)
                    acc = None; pivots0 = r['pivots']
                x = r['x']
                if acc is None:
                    acc = (list(x), r['p'])
                else:
                    xs, m = acc
                    xs = [core.crt_pair(a, m, b, r['p'])[0] for a, b in zip(xs, x)]
                    acc = (xs, m * r['p'])
                xs, m = acc
                rat = [core.ratrecon(a, m) for a in xs]
                if all(v is not None for v in rat):
                    N = [{} for _ in range(NV)]; Q = {}
                    for (kind, i, mon), v in zip(r['cols_list'], rat):
                        if v == 0: continue
                        if kind == 'N': N[i][mon] = v
                        elif kind == 'Q': Q[mon] = v
                    ok, info = core.verify_witness(F, A, q, D, N, Q)
                    if ok:
                        pool.terminate()
                        base.update({'N': N, 'Q': Q, 'primes_used': done})
                        return 'solved', base
                    elif verbose:
                        print(f'[parallel] reconstruction failed verification after {done} primes; continuing', file=sys.stderr)
    last = last or {}
    last['primes_used'] = done
    return 'lift-failed', last
