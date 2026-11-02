"""Multiprocessing sweep of candidate denominators for one instance.

Usage
-----
  python3 -m hodge.sweep --instance ID --candidates cands.jsonl --out results.jsonl \\
        [--workers 4] [--dmax 8] [--witness-dir DIR] [--log /tmp/hodge/results/ID/log.jsonl] \\
        [--timeout 1800] [--powers 1] [--no-e1] [--symmetric] [--dense]

Candidate file: one JSON object per line, {"D": "<expr>", "family": "<tag>"} (a bare
expression string per line is also accepted).  Each candidate is parsed with
core.parse; D must be homogeneous with pure Y-parity.

Per candidate (in a separate worker process, killed on timeout):
  1. degree cap and parity checks, F | D flag (jacobian.divisible_by),
  2. e1_solvable (A*D in J(F) mod p) unless --no-e1; if it fails the full system is skipped
     (e1 is a necessary condition, so the record says consistent=false),
  3. fixed-D solve: fixed_d.solve_fixed_D_fast (compressed rows + exact fmpq solve), with
     fallback to fixed_d.solve_fixed_D (dense rref + CRT lift) on exception or lift failure;
     --symmetric uses symmetry.solve_fixed_D_symmetric (character-sector restricted system)
     and --dense uses solve_fixed_D directly.
  4. On an exact witness, core.verify_witness is re-run and the witness is written to
     <witness-dir>/witness-<hash>.json in the witness.json format of the repo.

Output: one JSON line per candidate appended to --out and (if given) to --log with the
schema  {"D","d","family","kind","e1","consistent","verify","control","rank","nullity","solver","seconds","note", ...}
(kind: fixedD | e1 | skipped; verify: True only after core.verify_witness on an exact witness; control: from
--control or a "control": true candidate field).
Resumable: candidates whose normalized D already appears in --out are skipped.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys, time, traceback
import multiprocessing as mp
from . import core, fixed_d


def normalize_expr(expr):
    return core.to_str(core.parse(expr))


def read_candidates(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('{'):
                obj = json.loads(line)
            else:
                obj = {'D': line}
            obj.setdefault('family', 'unspecified')
            out.append(obj)
    return out


def read_done(path):
    done = set()
    if not path or not os.path.exists(path):
        return done
    with open(path) as f:
        for line in f:
            try:
                obj = json.loads(line)
            except Exception:
                continue
            if 'D_norm' in obj:
                done.add(obj['D_norm'])
    return done


def witness_record(iid, D, N, Q, family=None):
    return {'instance': iid, 'd': core.degree(D), 'D': core.to_str(D),
            'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q), 'engine': fixed_d.ENGINE_VERSION,
            'family': family}


def run_candidate(iid, F, A, q, cand, dmax, do_e1=True, mode='fast', symmetric=False, control=False):
    """Evaluate one candidate.  Returns the result record (dict).

    Record schema additions (engine v3): kind = 'fixedD' (full system solved) | 'e1' (stopped at the e1
    prefilter) | 'skipped' (parse/parity/degree-cap); verify = True only after core.verify_witness passed on
    an exact witness; control = True for calibration runs (candidate field "control": true or --control)."""
    t0 = time.time()
    rec = {'D': cand['D'], 'family': cand.get('family', 'unspecified'), 'instance': iid,
           'engine': fixed_d.ENGINE_VERSION, 'kind': 'fixedD', 'verify': False,
           'control': bool(control or cand.get('control', False)),
           'e1': None, 'consistent': None, 'rank': None, 'nullity': None, 'note': ''}
    try:
        D = core.parse(cand['D'])
    except Exception as ex:
        rec.update({'d': None, 'note': f'parse error: {ex}', 'seconds': round(time.time() - t0, 2)})
        return rec
    rec['D_norm'] = core.to_str(D)
    d = core.degree(D)
    rec['d'] = d
    if not D or not core.is_homogeneous(D):
        rec['note'] = 'D zero or not homogeneous'
    elif core.parity_of(D) is None:
        rec['note'] = 'D has mixed Y-parity'
    elif dmax is not None and d > dmax:
        rec['note'] = f'skipped: degree {d} > cap {dmax}'
    if rec['note']:
        rec['kind'] = 'skipped'
        rec['seconds'] = round(time.time() - t0, 2)
        return rec
    notes = []
    try:
        from . import jacobian
        if d >= 6:
            divF, _ = jacobian.divisible_by(F, D)
            if divF:
                notes.append('F divides D (trivial/spurious denominator)')
                rec['F_divides_D'] = True
    except Exception as ex:
        notes.append(f'divisibility check failed: {ex}')
    if do_e1:
        try:
            ok, info = fixed_d.e1_solvable(F, A, q, D)
            rec['e1'] = bool(ok)
            rec['e1_info'] = {k: info[k] for k in ('rank', 'rank_aug', 'rows', 'cols')}
            if not ok:
                rec['consistent'] = False
                rec['kind'] = 'e1'
                notes.append('e1 fails: A*D not in J(F) mod p')
                rec['note'] = '; '.join(notes)
                rec['seconds'] = round(time.time() - t0, 2)
                return rec
        except Exception as ex:
            notes.append(f'e1 error: {ex}')
    exact = None
    try:
        if symmetric:
            from . import symmetry
            status, data = symmetry.solve_fixed_D_symmetric(F, A, q, D)
            rec.update({'consistent': status != 'inconsistent', 'rank': data.get('rank'),
                        'nullity': data.get('nullity'), 'rows': data.get('rows'), 'cols': data.get('cols'),
                        'full_cols': data.get('full_cols'), 'subgroup_size': data.get('subgroup_size'),
                        'solver': 'symmetric'})
            if status == 'solved':
                exact = (data['N'], data['Q'])
            elif status == 'lift-failed':
                notes.append('symmetric lift failed')
        elif mode == 'dense':
            status, data = fixed_d.solve_fixed_D(F, A, q, D)
            rec.update({'consistent': status != 'inconsistent', 'rank': data.get('rank'),
                        'nullity': data.get('nullity'), 'rows': data.get('rows'), 'cols': data.get('cols'),
                        'solver': 'dense'})
            if status == 'solved':
                exact = (data['N'], data['Q'])
            elif status == 'lift-failed':
                notes.append('dense lift failed')
        else:
            try:
                out = fixed_d.solve_fixed_D_fast(F, A, q, D)
                rec.update({'consistent': bool(out['consistent']), 'rank': out['rank'],
                            'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
                            'solver': 'fast:' + str(out.get('solver', 'compressed'))})
                if out.get('sparse_stats'):
                    rec['sparse_stats'] = out['sparse_stats']
                if out['consistent']:
                    if out.get('exact'):
                        exact = (out['N'], out['Q'])
                    else:
                        notes.append('fast exact solve failed; falling back to dense CRT lift')
                        status, data = fixed_d.solve_fixed_D(F, A, q, D)
                        rec['solver'] = 'fast+dense'
                        if status == 'solved':
                            exact = (data['N'], data['Q'])
                        elif status == 'inconsistent':
                            rec['consistent'] = False
                            notes.append('dense solver says inconsistent (compressed test was wrong?)')
                        else:
                            notes.append('dense lift failed too')
            except Exception as ex:
                notes.append(f'fast solver error ({ex.__class__.__name__}: {ex}); fallback to dense')
                status, data = fixed_d.solve_fixed_D(F, A, q, D)
                rec.update({'consistent': status != 'inconsistent', 'rank': data.get('rank'),
                            'nullity': data.get('nullity'), 'rows': data.get('rows'),
                            'cols': data.get('cols'), 'solver': 'dense-fallback'})
                if status == 'solved':
                    exact = (data['N'], data['Q'])
    except Exception as ex:
        notes.append(f'solver error: {ex.__class__.__name__}: {ex}')
        rec['traceback'] = traceback.format_exc()[-2000:]
    if exact is not None:
        N, Q = exact
        ok, info = core.verify_witness(F, A, q, D, N, Q)
        rec['exact'] = bool(ok)
        rec['verify'] = bool(ok)
        if ok:
            rec['witness'] = witness_record(iid, D, N, Q, rec['family'])
            rec['strict'] = not Q
            notes.append('EXACT WITNESS verified' + (' (strict, Q=0)' if not Q else ''))
        else:
            notes.append(f'exact solution FAILED verification: {info}')
    rec['note'] = '; '.join(notes)
    rec['seconds'] = round(time.time() - t0, 2)
    return rec


# ----------------------------------------------------------------------------
# worker processes
# ----------------------------------------------------------------------------

def _worker(iid, cand, dmax, do_e1, mode, symmetric, conn, control=False):
    try:
        F, A, q, _ = core.load_instance(iid)
        rec = run_candidate(iid, F, A, q, cand, dmax, do_e1, mode, symmetric, control=control)
    except Exception as ex:
        rec = {'D': cand.get('D'), 'family': cand.get('family'), 'instance': iid, 'd': None,
               'e1': None, 'consistent': None, 'rank': None, 'nullity': None, 'seconds': None,
               'note': f'worker crashed: {ex.__class__.__name__}: {ex}'}
    try:
        conn.send(rec)
    except Exception:
        pass
    conn.close()


def sweep(iid, cands, out_path, workers=2, dmax=None, log_path=None, witness_dir=None,
          timeout=None, do_e1=True, mode='fast', symmetric=False, powers=(1,), verbose=True, control=False):
    """Run the sweep.  Returns list of result records (new ones only)."""
    done = read_done(out_path)
    jobs = []
    for c in cands:
        for k in powers:
            cc = dict(c)
            if k != 1:
                cc['D'] = f'({c["D"]})**{k}'
                cc['family'] = f'{c.get("family", "unspecified")}^{k}'
            try:
                key = normalize_expr(cc['D'])
            except Exception:
                key = cc['D']
            if key in done:
                continue
            done.add(key)
            jobs.append(cc)
    if verbose:
        print(f'[sweep] {iid}: {len(jobs)} candidates to run ({len(cands)} listed, powers {list(powers)}), '
              f'{workers} workers, dmax={dmax}, timeout={timeout}', file=sys.stderr)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    if log_path:
        os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
    if witness_dir:
        os.makedirs(witness_dir, exist_ok=True)
    ctx = mp.get_context('fork')
    running = []  # (proc, conn, cand, t_start)
    results = []
    idx = 0
    out_f = open(out_path, 'a')
    log_f = open(log_path, 'a') if log_path else None

    def emit(rec):
        results.append(rec)
        line = json.dumps(rec, default=str)
        out_f.write(line + '\n'); out_f.flush()
        if log_f:
            slim = {k: rec.get(k) for k in ('D', 'd', 'family', 'kind', 'e1', 'consistent', 'verify', 'control',
                                             'rank', 'nullity', 'solver', 'seconds', 'note')}
            slim['engine'] = rec.get('engine', fixed_d.ENGINE_VERSION)
            slim['instance'] = iid
            log_f.write(json.dumps(slim, default=str) + '\n'); log_f.flush()
        if rec.get('witness') and witness_dir:
            h = hashlib.sha1(rec['witness']['D'].encode()).hexdigest()[:10]
            wp = os.path.join(witness_dir, f'witness-{h}.json')
            with open(wp, 'w') as wf:
                json.dump(rec['witness'], wf, indent=2)
            if verbose:
                print(f'[sweep] *** EXACT WITNESS written to {wp}', file=sys.stderr)
        if verbose:
            print(f"[sweep] d={rec.get('d')} e1={rec.get('e1')} consistent={rec.get('consistent')} "
                  f"rank={rec.get('rank')} nullity={rec.get('nullity')} {rec.get('seconds')}s "
                  f"{rec.get('family')} :: {str(rec.get('D'))[:70]} :: {rec.get('note')}", file=sys.stderr)

    try:
        while idx < len(jobs) or running:
            while idx < len(jobs) and len(running) < workers:
                cand = jobs[idx]; idx += 1
                parent, child = ctx.Pipe(duplex=False)
                pr = ctx.Process(target=_worker, args=(iid, cand, dmax, do_e1, mode, symmetric, child, control))
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
                    emit({'D': cand['D'], 'family': cand.get('family'), 'instance': iid, 'd': None,
                          'e1': None, 'consistent': None, 'rank': None, 'nullity': None,
                          'seconds': round(time.time() - ts, 2), 'note': f'worker died (exit {pr.exitcode})'})
                    continue
                if timeout and time.time() - ts > timeout:
                    pr.terminate(); pr.join(5); conn.close()
                    emit({'D': cand['D'], 'family': cand.get('family'), 'instance': iid, 'd': None,
                          'e1': None, 'consistent': None, 'rank': None, 'nullity': None,
                          'seconds': round(time.time() - ts, 2), 'note': f'timeout after {timeout}s (killed)'})
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
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--instance', required=True)
    ap.add_argument('--candidates', required=True, help='JSONL of {"D": expr, "family": tag} or bare expressions')
    ap.add_argument('--out', required=True, help='results JSONL (appended; used for resume)')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--dmax', type=int, default=None, help='skip candidates of degree above this cap')
    ap.add_argument('--log', default=None, help='log.jsonl in the team schema (appended)')
    ap.add_argument('--witness-dir', default=None, help='directory for exact witnesses found')
    ap.add_argument('--timeout', type=float, default=None, help='seconds per candidate before the worker is killed')
    ap.add_argument('--powers', default='1', help='comma list of powers of each candidate to test, e.g. 1,2 at q=2')
    ap.add_argument('--no-e1', action='store_true', help='skip the e1 prefilter')
    ap.add_argument('--dense', action='store_true', help='use the dense rref + CRT solver instead of the compressed one')
    ap.add_argument('--symmetric', action='store_true', help='use the diagonal-symmetry sector-restricted solver')
    ap.add_argument('--control', action='store_true', help='mark every record control=true (calibration / illegal-D runs)')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    cands = read_candidates(a.candidates)
    powers = tuple(int(x) for x in a.powers.split(',') if x.strip())
    res = sweep(a.instance, cands, a.out, workers=a.workers, dmax=a.dmax, log_path=a.log,
                witness_dir=a.witness_dir, timeout=a.timeout, do_e1=not a.no_e1,
                mode='dense' if a.dense else 'fast', symmetric=a.symmetric, powers=powers,
                verbose=not a.quiet, control=a.control)
    nwit = sum(1 for r in res if r.get('witness'))
    ncons = sum(1 for r in res if r.get('consistent'))
    print(f'[sweep] done: {len(res)} run, {ncons} consistent, {nwit} exact witnesses', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
