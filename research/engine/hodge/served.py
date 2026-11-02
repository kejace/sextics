"""Served-space tool: which Y-odd targets A' of degree 6q does a denominator D serve?  (task 3)

For fixed (F, q, D) the cleared identity is linear in (N, Q, A'):  M_NQ (n, q) + M_A a = 0 with the A'
columns m * D^2.  The served space  T(D) = {a : exists (n,q)}  is the projection of the nullspace onto the
A' coordinates.  Here it is computed as the nullspace of the Schur constraints after eliminating the (N, Q)
unknowns sparsely (sge.schur_constraints), which avoids the full nullspace of fixed_d.target_image.

    from hodge import served
    r = served.served_space(F, q, D)        # dict(dim, basis, constraints, contains(A) helper, stats)
    served.contains(r, A)                   # is the instance target served by D?

CLI (one pass over a candidate list, multiprocessing, resumable, JSONL out):
    python3 -m hodge.served --instance ID --candidates cands.jsonl --out served.jsonl [--workers 4]
        [--dmax 10] [--timeout 3600] [--basis] [--log /tmp/hodge/results/ID/log.jsonl]
Candidate lines: {"D": "<expr>", "family": "<tag>"} or bare expressions.  Records:
    {"D","D_norm","d","family","kind":"served","dim","contains_A","nA","rank_constraints","seconds","note",
     "basis": [...] (with --basis; polys mod p in to_str form), "engine", "p"}
`dim` is the dimension of the served space mod p (an upper bound for characteristic 0; equality for all but
finitely many primes).  contains_A == True  <=>  the fixed-D system for the instance target A is consistent mod p.
"""
from __future__ import annotations
import argparse, json, os, sys, time, traceback
import multiprocessing as mp
import numpy as np
from . import core, fixed_d, sge
from .core import NV, monomials


def served_space(F, q, D, p=None, verbose=False, want_basis=True, **kw):
    """Served space of D mod p.  Returns dict(dim, nA, basis (list of polys mod p) or None, constraints (R),
    a_mons (monomial list of the A' coordinates), p, stats, S)."""
    p = p or core.PRIMES[0]
    A0 = {monomials(6 * q, 1)[0]: 1}
    S = fixed_d.FixedDSystem(F, A0, q, D, p, with_target_space=True, verbose=verbose)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    a_mons = [S.cols[j][2] for j in a_idx]
    R, stats = sge.schur_constraints(M, p, a_idx, verbose=verbose, **kw)
    nA = len(a_idx)
    rank_c = R.shape[0]
    out = {'dim': nA - rank_c, 'nA': nA, 'rank_constraints': rank_c, 'constraints': R, 'a_mons': a_mons, 'p': p,
           'stats': stats, 'S': S, 'basis': None}
    if want_basis and nA - rank_c > 0:
        from flint import nmod_mat
        if rank_c == 0:
            basis = [{m: 1} for m in a_mons]
        else:
            Rn = nmod_mat(rank_c, nA, R.ravel().tolist(), p)
            X, nul = Rn.nullspace()
            basis = []
            for k in range(nul):
                poly = {a_mons[j]: int(X[j, k]) for j in range(nA) if int(X[j, k])}
                basis.append(poly)
        out['basis'] = basis
    return out


def contains(res, A):
    """Is the target A (exact or mod p poly, Y-odd degree 6q) in the served space?  (R a == 0 mod p)"""
    p = res['p']
    Ap = core.reduce_mod(A, p)
    idx = {m: j for j, m in enumerate(res['a_mons'])}
    a = np.zeros(len(res['a_mons']), dtype=np.int64)
    for m, v in Ap.items():
        if m not in idx:
            return False
        a[idx[m]] = v
    R = res['constraints']
    if R.shape[0] == 0:
        return True
    lo, hi = a & 0xFFFF, a >> 16
    r = ((R @ lo) % p + ((R @ hi) % p) * (1 << 16)) % p
    return bool(np.all(r == 0))


def run_candidate(iid, F, A, q, cand, dmax=None, want_basis=False):
    t0 = time.time()
    rec = {'D': cand['D'], 'family': cand.get('family', 'unspecified'), 'instance': iid, 'kind': 'served',
           'engine': fixed_d.ENGINE_VERSION, 'dim': None, 'contains_A': None, 'note': ''}
    try:
        D = core.parse(cand['D'])
    except Exception as ex:
        rec.update({'d': None, 'note': f'parse error: {ex}', 'seconds': round(time.time() - t0, 2)})
        return rec
    rec['D_norm'] = core.to_str(D)
    d = core.degree(D); rec['d'] = d
    if not D or not core.is_homogeneous(D):
        rec['note'] = 'D zero or not homogeneous'
    elif core.parity_of(D) is None:
        rec['note'] = 'D has mixed Y-parity'
    elif dmax is not None and d > dmax:
        rec['note'] = f'skipped: degree {d} > cap {dmax}'
    if rec['note']:
        rec['seconds'] = round(time.time() - t0, 2)
        return rec
    try:
        r = served_space(F, q, D, want_basis=want_basis)
        rec.update({'dim': r['dim'], 'nA': r['nA'], 'rank_constraints': r['rank_constraints'], 'p': r['p'],
                    'contains_A': contains(r, A), 'rows': r['S'].nrows, 'cols': r['S'].ncols,
                    'sparse_stats': {k: r['stats'].get(k) for k in ('status_sparse', 'npiv_sparse', 'active_rows',
                                                                   'active_unknown_cols', 't_sparse', 't_total')}})
        if want_basis and r['basis'] is not None:
            rec['basis'] = [core.to_str({m: (v if v <= r['p'] // 2 else v - r['p']) for m, v in b.items()}) for b in r['basis']]
    except Exception as ex:
        rec['note'] = f'error: {ex.__class__.__name__}: {ex}'
        rec['traceback'] = traceback.format_exc()[-1500:]
    rec['seconds'] = round(time.time() - t0, 2)
    return rec


def _worker(iid, cand, dmax, want_basis, conn):
    os.environ['HODGE_QUIET'] = '1'
    try:
        F, A, q, _ = core.load_instance(iid)
        rec = run_candidate(iid, F, A, q, cand, dmax, want_basis)
    except Exception as ex:
        rec = {'D': cand.get('D'), 'family': cand.get('family'), 'instance': iid, 'kind': 'served', 'd': None,
               'dim': None, 'contains_A': None, 'seconds': None, 'note': f'worker crashed: {ex.__class__.__name__}: {ex}'}
    try:
        conn.send(rec)
    except Exception:
        pass
    conn.close()


def sweep(iid, cands, out_path, workers=2, dmax=None, log_path=None, timeout=None, want_basis=False, verbose=True):
    from .sweep import read_done, normalize_expr
    done = read_done(out_path)
    jobs = []
    for c in cands:
        try:
            key = normalize_expr(c['D'])
        except Exception:
            key = c['D']
        if key in done:
            continue
        done.add(key); jobs.append(c)
    if verbose:
        print(f'[served] {iid}: {len(jobs)} candidates to run ({len(cands)} listed), {workers} workers', file=sys.stderr)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    ctx = mp.get_context('fork')
    running = []; results = []; idx = 0
    out_f = open(out_path, 'a'); log_f = open(log_path, 'a') if log_path else None

    def emit(rec):
        results.append(rec)
        out_f.write(json.dumps(rec, default=str) + '\n'); out_f.flush()
        if log_f:
            slim = {k: rec.get(k) for k in ('D', 'd', 'family', 'kind', 'dim', 'nA', 'contains_A', 'seconds', 'note')}
            slim.update({'engine': fixed_d.ENGINE_VERSION, 'instance': iid, 'consistent': rec.get('contains_A'),
                         'verify': False})
            log_f.write(json.dumps(slim, default=str) + '\n'); log_f.flush()
        if verbose:
            print(f"[served] d={rec.get('d')} dim={rec.get('dim')}/{rec.get('nA')} contains_A={rec.get('contains_A')} "
                  f"{rec.get('seconds')}s {rec.get('family')} :: {str(rec.get('D'))[:70]} :: {rec.get('note')}", file=sys.stderr)
    try:
        while idx < len(jobs) or running:
            while idx < len(jobs) and len(running) < workers:
                cand = jobs[idx]; idx += 1
                parent, child = ctx.Pipe(duplex=False)
                pr = ctx.Process(target=_worker, args=(iid, cand, dmax, want_basis, child))
                pr.start(); child.close()
                running.append((pr, parent, cand, time.time()))
            still = []
            for pr, conn, cand, ts in running:
                got = None
                if conn.poll(0.2):
                    try:
                        got = conn.recv()
                    except EOFError:
                        got = None
                if got is not None:
                    pr.join(); conn.close(); emit(got); continue
                if not pr.is_alive():
                    pr.join(); conn.close()
                    emit({'D': cand['D'], 'family': cand.get('family'), 'instance': iid, 'kind': 'served', 'd': None,
                          'dim': None, 'contains_A': None, 'seconds': round(time.time() - ts, 2),
                          'note': f'worker died (exit {pr.exitcode})'})
                    continue
                if timeout and time.time() - ts > timeout:
                    pr.terminate(); pr.join(5); conn.close()
                    emit({'D': cand['D'], 'family': cand.get('family'), 'instance': iid, 'kind': 'served', 'd': None,
                          'dim': None, 'contains_A': None, 'seconds': round(time.time() - ts, 2),
                          'note': f'timeout after {timeout}s (killed)'})
                    continue
                still.append((pr, conn, cand, ts))
            running = still
    finally:
        for pr, conn, cand, ts in running:
            pr.terminate()
        out_f.close()
        if log_f:
            log_f.close()
    return results


def main(argv=None):
    from .sweep import read_candidates
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--instance', required=True)
    ap.add_argument('--candidates', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--dmax', type=int, default=None)
    ap.add_argument('--log', default=None)
    ap.add_argument('--timeout', type=float, default=None)
    ap.add_argument('--basis', action='store_true', help='include a basis of the served space (mod p) in each record')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    res = sweep(a.instance, read_candidates(a.candidates), a.out, workers=a.workers, dmax=a.dmax, log_path=a.log,
                timeout=a.timeout, want_basis=a.basis, verbose=not a.quiet)
    n_in = sum(1 for r in res if r.get('contains_A'))
    print(f'[served] done: {len(res)} run, {n_in} serve the instance target', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
