"""Wave-2 E1-kernel re-verification sweep (independent prime + two code paths).

Usage: python3 w2_e1kernel.py <prime-index> <jobs...>     e.g.  python3 w2_e1kernel.py 0 1e 1o 2e ... 10o
Runs the tasks in a multiprocessing pool (one process per task).  Each task:
  * w2lib.e1_kernel (own implementation: rank of J_{d+6} and [J | A*S_d], rank of J_d) with the
    Hilbert-series check (R1);
  * hodge.jacobian.mult_kernel (engine code path) for the same (d, parity, prime);
  * appends one JSON line to wave2-excl.log.jsonl and to w2_e1kernel.<prime>.jsonl.
"""
import sys, time, json, os
import multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import w2lib as W
from hodge import jacobian

OUT = None


def task(args):
    d, par, p = args
    F, A, q, rec = W.load()
    t0 = time.time()
    info, _ = W.e1_kernel(F, A, d, par, p, want_basis=False)
    t1 = time.time()
    try:
        mk = jacobian.mult_kernel(F, A, d, p=p, parity=par)
        eng = {'engine_dim': mk['dim'], 'engine_kernel_dim': mk['kernel_dim'], 'engine_jacobian_dim': mk['jacobian_dim'],
               'engine_ring_dim': mk['ring_dim'], 'engine_seconds': round(time.time() - t1, 1)}
    except Exception as ex:
        eng = {'engine_error': repr(ex)}
    rec = {'exp': 'e1-kernel', 'd': d, 'parD': par, 'prime': p, 'result': info['kernel mod J'],
           'rank': info['rank J_{d+6}'], 'rank_expected': info['rank J_{d+6} expected'], 'dimJ_d': info['dim J_d'],
           'dimJ_d_expected': info['dim J_d expected'], 'hilbert_ok': info['hilbert_ok'],
           'kernel_incl_J': info['kernel dim (incl J_d)'], 'rows': info['rows'], 'cols': info['cols'],
           'seconds': round(time.time() - t0, 1), 'job': os.environ.get('W2JOB', 'w2_e1kernel'),
           'note': 'dim_Q(K_d/J_d) <= result if hilbert_ok (R1+R2); engine cross-check: ' + json.dumps(eng)}
    rec.update(eng)
    W.wlog(**rec)
    with open(OUT, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec


if __name__ == '__main__':
    pidx = int(sys.argv[1])
    p = W.PRIMES30[pidx] if pidx >= 0 else W.PRIMES31[-pidx - 1]
    OUT = os.path.join(W.RES, f'w2_e1kernel.p{p}.jsonl')
    jobs = []
    for j in sys.argv[2:]:
        d = int(j[:-1]); par = 0 if j[-1] == 'e' else 1
        jobs.append((d, par, p))
    nproc = min(len(jobs), int(os.environ.get('W2NPROC', '16')))
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(task, jobs):
            print('done', r['d'], r['parD'], r['result'], flush=True)
