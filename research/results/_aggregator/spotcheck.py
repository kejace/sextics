"""Aggregator spot-checks: independent recomputation of E1-kernel dimensions (Ann_R(A)_d) at a fresh prime.
One JSONL line per (lane, d, parity) appended to spotcheck.jsonl as soon as it finishes."""
import sys, json, time, os, traceback
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from multiprocessing import Pool
from hodge import core, jacobian
P = 2147483563  # core.PRIMES[4]; not used by any lane so far
OUT = '/tmp/hodge/results/_aggregator/spotcheck.jsonl'
# (lane, d, parity, expected-from-reports)
TASKS = [
 ('symmetric-sextic-residual',6,0,'0'),('symmetric-sextic-residual',6,1,'0'),('symmetric-sextic-residual',7,0,'1 (D7)'),('symmetric-sextic-residual',7,1,'1 (D7Y)'),
 ('short-power-presentation',6,0,'0'),('short-power-presentation',6,1,'0'),('short-power-presentation',7,0,'1 (wave-1)'),('short-power-presentation',7,1,'1 (wave-1)'),
 ('rank-one-biquartic',6,0,'0'),('rank-one-biquartic',6,1,'0'),('rank-one-biquartic',7,0,'3'),('rank-one-biquartic',7,1,'0'),
 ('one-reflection-skew',6,0,'0'),('one-reflection-skew',6,1,'0'),('one-reflection-skew',7,0,'? (total 9)'),('one-reflection-skew',7,1,'? (total 9)'),
 ('off-paired',4,0,'1 (D0)'),('off-paired',4,1,'0'),('off-paired',5,0,'? (total 5)'),('off-paired',5,1,'? (total 5)'),('off-paired',6,0,'? (total 18)'),('off-paired',6,1,'? (total 18)'),('off-paired',7,0,'? (total 51)'),('off-paired',7,1,'? (total 51)'),
 ('four-plus-two-residual',5,0,'0'),('four-plus-two-residual',5,1,'0'),('four-plus-two-residual',6,0,'20'),('four-plus-two-residual',6,1,'13 (core odd)'),('four-plus-two-residual',7,0,'63'),
 ('cyclic-triple-residual',7,0,'0'),('cyclic-triple-residual',7,1,'0'),('cyclic-triple-residual',8,0,'1 (D8)'),('cyclic-triple-residual',8,1,'0'),
 ('diagonal-permutation',7,0,'0'),('diagonal-permutation',7,1,'0'),('diagonal-permutation',8,0,'1 (D8)'),('diagonal-permutation',8,1,'0'),
 ('quadric-contact-residual',7,0,'0'),('quadric-contact-residual',7,1,'0'),('quadric-contact-residual',8,0,'? (total 1: D8*)'),('quadric-contact-residual',8,1,'? (total 1)'),
 ('reflection-coupled-residual',7,0,'0'),('reflection-coupled-residual',7,1,'0'),('reflection-coupled-residual',8,0,'1 (K8)'),('reflection-coupled-residual',8,1,'? (not computed)'),
 ('quadratic-cubic-pullback',7,0,'0'),('quadratic-cubic-pullback',7,1,'0'),('quadratic-cubic-pullback',8,0,'0'),('quadratic-cubic-pullback',8,1,'? (not computed)'),
 ('klein-coupled-character',8,0,'0'),('klein-coupled-character',8,1,'? (not computed)'),
 ('polarized-k3-remainder',8,0,'0'),('polarized-k3-remainder',8,1,'0'),
 ('dense-rational',8,0,'0'),('dense-rational',8,1,'0'),
 ('cycle-q2',4,0,'0'),('cycle-q2',4,1,'0'),('cycle-q2',5,0,'0'),('cycle-q2',5,1,'0'),
 ('symmetric-sextic-residual',8,0,'6'),('symmetric-sextic-residual',8,1,'6'),
 ('short-power-presentation',8,0,'6 (wave-1)'),('short-power-presentation',8,1,'? (not computed)'),
 ('rank-one-biquartic',8,0,'9'),('rank-one-biquartic',8,1,'9'),
 ('one-reflection-skew',8,0,'? (total 144)'),('one-reflection-skew',8,1,'? (total 144)'),
 ('off-paired',8,0,'? (total 121)'),('off-paired',8,1,'? (total 121)'),
 ('eight-plane-configuration',9,0,'0'),('eight-plane-configuration',9,1,'0'),
 ('dense-rational',9,0,'0'),('dense-rational',9,1,'0'),
 ('polarized-k3-remainder',9,0,'0'),('polarized-k3-remainder',9,1,'0'),
 ('klein-coupled-character',9,0,'24'),('klein-coupled-character',9,1,'24'),
 ('diagonal-permutation',9,0,'? (unknown k_9)'),('diagonal-permutation',9,1,'? (unknown k_9)'),
 ('quadric-contact-residual',9,0,'? (total 4)'),('quadric-contact-residual',9,1,'? (total 4)'),
 ('quadratic-cubic-pullback',9,0,'? (running elsewhere)'),('quadratic-cubic-pullback',9,1,'?'),
 ('cyclic-triple-residual',9,0,'? (total 120)'),('cyclic-triple-residual',9,1,'? (total 120)'),
 ('reflection-coupled-residual',9,0,'? (>=21 forced (+,+))'),('reflection-coupled-residual',9,1,'?'),
 ('four-plus-two-residual',7,1,'?'),('four-plus-two-residual',8,0,'?'),('four-plus-two-residual',8,1,'?'),
 ('rank-one-biquartic',9,0,'27'),('rank-one-biquartic',9,1,'27'),
 ('symmetric-sextic-residual',9,0,'?'),('symmetric-sextic-residual',9,1,'?'),
 ('short-power-presentation',9,0,'?'),('short-power-presentation',9,1,'?'),
 ('one-reflection-skew',9,0,'? (total 684)'),('one-reflection-skew',9,1,'? (total 684)'),
]
_cache = {}
def run(task):
    lane, d, par, exp = task
    t0 = time.time()
    try:
        F, A, q, rec = core.load_instance(lane, P)
        r = jacobian.mult_kernel(F, A, d, P, parity=par)
        rec_out = {'lane': lane, 'd': d, 'parity': par, 'p': P, 'q': q, 'degA': core.degree(A),
                   'ann_mod_J': r['dim'], 'kernel_dim': r['kernel_dim'], 'jacobian_dim': r['jacobian_dim'],
                   'ring_dim': r['ring_dim'], 'rows': r['rows'], 'cols': r['cols'],
                   'expected': exp, 'seconds': round(time.time()-t0,1), 'engine': 2,
                   'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        # consistency: implied rank of J block = cols - nullity ... report nullity-derived J rank
        rec_out['implied_rank_J_target'] = r['cols'] - len(r['basis']) if False else None
    except Exception as e:
        rec_out = {'lane': lane, 'd': d, 'parity': par, 'p': P, 'error': repr(e), 'trace': traceback.format_exc()[-800:],
                   'seconds': round(time.time()-t0,1), 'ts': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    with open(OUT, 'a') as f:
        f.write(json.dumps(rec_out) + '\n')
    return rec_out
if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                r = json.loads(line)
                if 'ann_mod_J' in r: done.add((r['lane'], r['d'], r['parity']))
            except Exception: pass
    todo = [t for t in TASKS if (t[0], t[1], t[2]) not in done]
    print(f'{len(todo)} tasks, {nproc} workers', flush=True)
    with Pool(nproc) as pool:
        for r in pool.imap_unordered(run, todo):
            print(json.dumps({k: r.get(k) for k in ('lane','d','parity','ann_mod_J','expected','seconds','error')}), flush=True)
