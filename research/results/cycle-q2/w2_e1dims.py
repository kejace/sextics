"""Step 1: E1-admissible denominator spaces per degree and character sector (mod p, two primes on request),
K_kappa = {D in (S_d)_kappa : A D in J(F)}, with dim (J_d)_kappa and dim (S_d)_kappa.
Usage: python3 w2_e1dims.py DMIN DMAX [NPROC] [PRIME]
Output: w2_e1dims_d{d}.json per degree (dims + mod-p kernel bases), log lines in wave2-search.log.jsonl."""
import os as _os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "FLINT_NUM_THREADS"):
    _os.environ[_v] = "1"
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import multiprocessing as mp
import w2lib as W
from hodge import core

OUT = '/tmp/hodge/results/cycle-q2'


def work(args):
    d, kappa, p = args
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    t0 = time.time()
    r = W.e1_kernel_modp(S, d, kappa, q, p=p)
    par = S.parity_of_key(d, kappa)
    basis = [core.to_str({e: v for e, v in b.items()}) for b in r['basis']]
    return {'d': d, 'key': kappa, 'parity': par, 'dimS': r['dimS'], 'dimK': r['dimK'], 'dimJ': r['dimJ'],
            'dimK_mod_J': r['dimK'] - r['dimJ'], 'p': p, 'seconds': round(time.time() - t0, 2), 'basis_modp': basis}


if __name__ == '__main__':
    dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    p = int(sys.argv[4]) if len(sys.argv) > 4 else W.PRIME
    F, A, q, rec = core.load_instance('cycle-q2')
    S = W.Sectors(F, A)
    for d in range(dmin, dmax + 1):
        t0 = time.time()
        keys = S.all_keys(d)
        with mp.get_context('fork').Pool(nproc) as pool:
            res = pool.map(work, [(d, k, p) for k in keys])
        res.sort(key=lambda r: r['key'])
        tot = {'d': d, 'p': p, 'n_keys': len(keys), 'dimS': sum(r['dimS'] for r in res), 'dimK': sum(r['dimK'] for r in res),
               'dimJ': sum(r['dimJ'] for r in res), 'dimK_mod_J': sum(r['dimK_mod_J'] for r in res),
               'dimS_even': sum(r['dimS'] for r in res if r['parity'] == 0), 'dimK_even': sum(r['dimK'] for r in res if r['parity'] == 0),
               'dimJ_even': sum(r['dimJ'] for r in res if r['parity'] == 0),
               'per_key': {str(r['key']): [r['parity'], r['dimS'], r['dimK'], r['dimJ']] for r in res},
               'seconds': round(time.time() - t0, 1)}
        with open(f'{OUT}/w2_e1dims_d{d}.json', 'w') as f:
            json.dump({'summary': tot, 'sectors': res}, f)
        W.log(step='e1dims', d=d, p=p, dimS=tot['dimS'], dimK=tot['dimK'], dimJ=tot['dimJ'], dimK_mod_J=tot['dimK_mod_J'],
              dimS_even=tot['dimS_even'], dimK_even=tot['dimK_even'], dimJ_even=tot['dimJ_even'],
              nonzero_sectors=sum(1 for r in res if r['dimK'] > r['dimJ']), seconds=tot['seconds'],
              note=f'E1 kernel per sector d={d}: total dimK={tot["dimK"]} (J {tot["dimJ"]}, new {tot["dimK_mod_J"]}) of dimS={tot["dimS"]}')
        print(f'd={d}: dimS={tot["dimS"]} dimK={tot["dimK"]} dimJ={tot["dimJ"]} new={tot["dimK_mod_J"]} '
              f'(even: S {tot["dimS_even"]} K {tot["dimK_even"]} J {tot["dimJ_even"]}) [{tot["seconds"]}s]', flush=True)
        big = sorted(res, key=lambda r: -(r['dimK'] - r['dimJ']))[:6]
        print('   largest new kernels:', [(r['key'], r['parity'], r['dimS'], r['dimK'], r['dimJ']) for r in big], flush=True)
