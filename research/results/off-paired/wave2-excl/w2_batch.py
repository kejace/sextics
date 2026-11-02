"""Multiprocess batch of rigorous fixed-D receipts (E1 with exact Koszul dimension; full restricted solve with the
explicit-kernel rank sandwich) on two primes.

usage: python3 w2_batch.py --cands cands.jsonl --out receipts.jsonl --job NAME [--workers N] [--primes p1,p2]
       [--dmax 14] [--timeout 3600] [--e1-only] [--no-second-prime]
Each candidate: {"D": expr over block names, "family": tag, "label": name}.  Resumable (labels already in --out skipped).
"""
import sys, os, json, time, argparse, traceback
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from hodge import core

def work(args):
    cand, primes, e1_only, job = args
    t0 = time.time()
    base = {'exp': 'fixedD-batch', 'job': job, 'label': cand.get('label', cand['D']), 'family': cand.get('family', '')}
    try:
        F, A, q = w2lib.load()
        D = w2lib.parse_expr(cand['D'])
        if not D or not core.is_homogeneous(D):
            return [dict(base, result='skip: not homogeneous', seconds=0)]
        if core.parity_of(D) is None:
            return [dict(base, result='skip: mixed Y-parity', d=core.degree(D), seconds=0)]
        if w2lib.char_of(D) is None:
            return [dict(base, result='skip: not an eigenvector', d=core.degree(D), seconds=0)]
        recs = w2lib.full_receipt(F, A, q, D, primes, cand.get('label', cand['D']), cand.get('family', ''), do_full=not e1_only)
        out = []
        for r in recs:
            r = dict(r); r['job'] = job; r['label'] = base['label']; r['D_expr'] = cand['D']
            r.pop('x', None); r.pop('S', None)
            out.append(r)
        return out
    except Exception as ex:
        return [dict(base, result='error: ' + repr(ex), trace=traceback.format_exc()[-800:], seconds=round(time.time() - t0, 1))]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cands', required=True); ap.add_argument('--out', required=True); ap.add_argument('--job', required=True)
    ap.add_argument('--workers', type=int, default=8); ap.add_argument('--primes', default=f'{w2lib.P1},{w2lib.P2}')
    ap.add_argument('--dmax', type=int, default=14); ap.add_argument('--e1-only', action='store_true')
    ap.add_argument('--timeout', type=int, default=7200)
    a = ap.parse_args()
    primes = [int(x) for x in a.primes.split(',')]
    cands = [json.loads(l) for l in open(a.cands) if l.strip()]
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            try: done.add(json.loads(l).get('label'))
            except Exception: pass
    todo = [c for c in cands if c.get('label', c['D']) not in done]
    print(f'[{a.job}] {len(todo)} candidates to run ({len(done)} done), workers {a.workers}, primes {primes}', flush=True)
    t0 = time.time()
    with mp.get_context('fork').Pool(a.workers, maxtasksperchild=4) as pool:
        for recs in pool.imap_unordered(work, [(c, primes, a.e1_only, a.job) for c in todo]):
            with open(a.out, 'a') as f:
                for r in recs:
                    f.write(json.dumps(r, default=str) + '\n')
            for r in recs:
                w2lib.log(r)
            print(json.dumps({k: recs[-1].get(k) for k in ('label', 'd', 'char', 'exp', 'prime', 'rank_M', 'rank_aug', 'explicit_kernel', 'rank_J', 'dimJ_exact', 'result', 'seconds')}, default=str), flush=True)
    print(f'[{a.job}] done in {time.time() - t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
