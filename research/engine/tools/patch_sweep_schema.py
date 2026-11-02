"""Idempotent: sweep records get kind / verify / control fields; --control flag; solver name from the fast path."""
import sys
P = '/home/kejace/devel/sextics-work/hodge/sweep.py'
s = open(P).read()
if "'kind': 'fixedD'" in s:
    print('sweep schema: already patched'); sys.exit(0)
old = '''def run_candidate(iid, F, A, q, cand, dmax, do_e1=True, mode='fast', symmetric=False):
    """Evaluate one candidate.  Returns the result record (dict)."""
    t0 = time.time()
    rec = {'D': cand['D'], 'family': cand.get('family', 'unspecified'), 'instance': iid,
           'engine': fixed_d.ENGINE_VERSION, 'e1': None, 'consistent': None, 'rank': None,
           'nullity': None, 'note': ''}'''
new = '''def run_candidate(iid, F, A, q, cand, dmax, do_e1=True, mode='fast', symmetric=False, control=False):
    """Evaluate one candidate.  Returns the result record (dict).

    Record schema additions (engine v3): kind = 'fixedD' (full system solved) | 'e1' (stopped at the e1
    prefilter) | 'skipped' (parse/parity/degree-cap); verify = True only after core.verify_witness passed on
    an exact witness; control = True for calibration runs (candidate field "control": true or --control)."""
    t0 = time.time()
    rec = {'D': cand['D'], 'family': cand.get('family', 'unspecified'), 'instance': iid,
           'engine': fixed_d.ENGINE_VERSION, 'kind': 'fixedD', 'verify': False,
           'control': bool(control or cand.get('control', False)),
           'e1': None, 'consistent': None, 'rank': None, 'nullity': None, 'note': ''}'''
assert old in s; s = s.replace(old, new)
old = '''    if rec['note']:
        rec['seconds'] = round(time.time() - t0, 2)
        return rec
    notes = []'''
new = '''    if rec['note']:
        rec['kind'] = 'skipped'
        rec['seconds'] = round(time.time() - t0, 2)
        return rec
    notes = []'''
assert old in s; s = s.replace(old, new)
old = '''            if not ok:
                rec['consistent'] = False
                notes.append('e1 fails: A*D not in J(F) mod p')'''
new = '''            if not ok:
                rec['consistent'] = False
                rec['kind'] = 'e1'
                notes.append('e1 fails: A*D not in J(F) mod p')'''
assert old in s; s = s.replace(old, new)
old = '''                out = fixed_d.solve_fixed_D_fast(F, A, q, D)
                rec.update({'consistent': bool(out['consistent']), 'rank': out['rank'],
                            'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
                            'solver': 'fast'})'''
new = '''                out = fixed_d.solve_fixed_D_fast(F, A, q, D)
                rec.update({'consistent': bool(out['consistent']), 'rank': out['rank'],
                            'nullity': out['nullity'], 'rows': out['rows'], 'cols': out['cols'],
                            'solver': 'fast:' + str(out.get('solver', 'compressed'))})
                if out.get('sparse_stats'):
                    rec['sparse_stats'] = out['sparse_stats']'''
assert old in s; s = s.replace(old, new)
old = '''        rec['exact'] = bool(ok)
        if ok:
            rec['witness'] = witness_record(iid, D, N, Q, rec['family'])'''
new = '''        rec['exact'] = bool(ok)
        rec['verify'] = bool(ok)
        if ok:
            rec['witness'] = witness_record(iid, D, N, Q, rec['family'])'''
assert old in s; s = s.replace(old, new)
old = '''def _worker(iid, cand, dmax, do_e1, mode, symmetric, conn):
    try:
        F, A, q, _ = core.load_instance(iid)
        rec = run_candidate(iid, F, A, q, cand, dmax, do_e1, mode, symmetric)'''
new = '''def _worker(iid, cand, dmax, do_e1, mode, symmetric, conn, control=False):
    try:
        F, A, q, _ = core.load_instance(iid)
        rec = run_candidate(iid, F, A, q, cand, dmax, do_e1, mode, symmetric, control=control)'''
assert old in s; s = s.replace(old, new)
old = '''def sweep(iid, cands, out_path, workers=2, dmax=None, log_path=None, witness_dir=None,
          timeout=None, do_e1=True, mode='fast', symmetric=False, powers=(1,), verbose=True):'''
new = '''def sweep(iid, cands, out_path, workers=2, dmax=None, log_path=None, witness_dir=None,
          timeout=None, do_e1=True, mode='fast', symmetric=False, powers=(1,), verbose=True, control=False):'''
assert old in s; s = s.replace(old, new)
old = '''                pr = ctx.Process(target=_worker, args=(iid, cand, dmax, do_e1, mode, symmetric, child))'''
new = '''                pr = ctx.Process(target=_worker, args=(iid, cand, dmax, do_e1, mode, symmetric, child, control))'''
assert old in s; s = s.replace(old, new)
old = '''            slim = {k: rec.get(k) for k in ('D', 'd', 'family', 'e1', 'consistent', 'rank', 'nullity', 'seconds', 'note')}'''
new = '''            slim = {k: rec.get(k) for k in ('D', 'd', 'family', 'kind', 'e1', 'consistent', 'verify', 'control',
                                             'rank', 'nullity', 'solver', 'seconds', 'note')}'''
assert old in s; s = s.replace(old, new)
old = '''    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)'''
new = '''    ap.add_argument('--control', action='store_true', help='mark every record control=true (calibration / illegal-D runs)')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)'''
assert old in s; s = s.replace(old, new)
old = '''                mode='dense' if a.dense else 'fast', symmetric=a.symmetric, powers=powers,
                verbose=not a.quiet)'''
new = '''                mode='dense' if a.dense else 'fast', symmetric=a.symmetric, powers=powers,
                verbose=not a.quiet, control=a.control)'''
assert old in s; s = s.replace(old, new)
s = s.replace('''Output: one JSON line per candidate appended to --out and (if given) to --log with the
schema  {"D","d","family","e1","consistent","rank","nullity","seconds","note", ...}.''',
'''Output: one JSON line per candidate appended to --out and (if given) to --log with the
schema  {"D","d","family","kind","e1","consistent","verify","control","rank","nullity","solver","seconds","note", ...}
(kind: fixedD | e1 | skipped; verify: True only after core.verify_witness on an exact witness; control: from
--control or a "control": true candidate field).''')
open(P, 'w').write(s); print('sweep schema: patched')
