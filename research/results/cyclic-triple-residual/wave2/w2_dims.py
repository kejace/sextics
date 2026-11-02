"""Exact E1-admissible spaces per (degree, sector): J_d^c independent generators + rationally reconstructed
kernel representatives.  Writes wave2/e1space_d{d}_c{a}{b}.json and logs one batch line per degree.
usage: python3 w2_dims.py dmin dmax [workers]
"""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2')
from w2lib import *
from multiprocessing import Pool

dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
inst = Inst()


def work(args):
    d, c = args
    fn = f'{W2}/e1space_d{d}_c{c[0]}{c[1]}.json'
    if os.path.exists(fn):
        return d, c, json.load(open(fn)), True
    t0 = time.time()
    meta, Jb, reps = inst.e1_space(d, c, nprimes=4)
    meta['seconds_total'] = round(time.time() - t0, 1)
    meta['J_indep'] = [[i, list(m)] for (i, m) in meta['J_indep']]
    with open(fn, 'w') as f:
        json.dump(meta, f)
    return d, c, meta, False


if __name__ == '__main__':
    for d in range(dmin, dmax + 1):
        t0 = time.time()
        tasks = [(d, c) for c in inst.chars]
        with Pool(workers) as pool:
            res = pool.map(work, tasks)
        tot = {'dim_S': 0, 'dimJ': 0, 'kerA': 0, 'kmodJ': 0}
        sectors = {}
        for d_, c, meta, cached in res:
            tot['dim_S'] += meta['dim_S']; tot['dimJ'] += meta['dimJ']; tot['kerA'] += meta['kerA_in_S']; tot['kmodJ'] += meta['ker_mod_J']
            sectors[f'{c[0]}{c[1]}'] = {'dim_S': meta['dim_S'], 'dimJ': meta['dimJ'], 'ker_mod_J': meta['ker_mod_J'],
                                       'family_dim': meta['dimJ'] + meta['ker_mod_J'], 'yparity': c[1] % 2}
        secs = round(time.time() - t0, 1)
        line = {'batch': 'e1-dims', 'd': d, 'totals': tot, 'sectors': sectors, 'seconds': secs,
                'note': 'exact E1-admissible family per G-character (order 18); kernel reps rationally reconstructed from 4 primes and verified mod a 5th'}
        log_batch(line)
        print(json.dumps(line), flush=True)
