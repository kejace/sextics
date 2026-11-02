"""E1 kernel K_d^{par} with multi-prime rational reconstruction (+ optional rational numerators).
Usage: python3 w2_kernel.py <d> <par> [--primes N] [--threads T] [--numerators] [--tag name]
Writes kernel_d{d}_p{par}.json and appends to wave2-search.log.jsonl."""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('par', type=int)
    ap.add_argument('--primes', type=int, default=3); ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--numerators', action='store_true'); ap.add_argument('--tag', default='')
    ap.add_argument('--maxprimes', type=int, default=6)
    a = ap.parse_args()
    set_threads(a.threads)
    F, A, q, rec = load()
    d, par = a.d, a.par
    out_path = os.path.join(RES, f'kernel_d{d}_p{par}.json')
    t0 = time.time()
    bases = {}; info = {}; rat = None; rinfo = {}
    primes = (PGB + core.PRIMES)[:a.maxprimes]
    for k, p in enumerate(primes):
        res = e1_kernel(F, A, d, par, p, verbose=True)
        bases[p] = res['basis_modJ']
        info[p] = {kk: res[kk] for kk in ('dim_S_d', 'dim_J_d', 'dim_R_d', 'dim_J_d6', 'dim_R_d6', 'kernel_dim_incl_J',
                                           'kernel_mod_J', 'nP', 'nT', 't_rows', 't_null', 't_total')}
        info[p]['t_jref_build'] = round(res['jref'].t_build, 1); info[p]['t_jref_rref'] = round(res['jref'].t_rref, 1)
        # save progress after each prime
        save = {'d': d, 'par': par, 'primes': list(bases), 'info': info,
                'basis_modp': {str(pp): [modp_poly_to_json(P) for P in bases[pp]] for pp in bases},
                'basis_rational': None, 'elapsed': round(time.time() - t0, 1)}
        with open(out_path, 'w') as f: json.dump(save, f)
        dims = {info[pp]['kernel_mod_J'] for pp in bases}
        if len(dims) != 1:
            print(f'[kernel] WARNING: kernel dims differ between primes: {dims}', file=sys.stderr, flush=True)
        if k + 1 >= 2 and len(dims) == 1 and bases[p]:
            rat, rinfo = reconstruct_polys(bases)
            print(f'[kernel] reconstruction with {k+1} primes: {"ok" if rat else "failed"} {rinfo}', file=sys.stderr, flush=True)
            if rat:
                save['basis_rational'] = [poly_to_json(P) for P in rat]
                save['reconstruction'] = rinfo
                with open(out_path, 'w') as f: json.dump(save, f)
                break
        if k + 1 >= a.primes and (not bases[p]):
            break
        if k + 1 >= a.primes and rat is None and k + 1 >= a.maxprimes:
            break
    kd = info[primes[0]]['kernel_mod_J']
    logrec = {'kind': 'e1-kernel', 'd': d, 'parity': par, 'kernel_mod_J': kd,
              'kernel_dim_incl_J': info[primes[0]]['kernel_dim_incl_J'], 'dim_J_d': info[primes[0]]['dim_J_d'],
              'dim_R_d': info[primes[0]]['dim_R_d'], 'dim_R_d6': info[primes[0]]['dim_R_d6'],
              'primes': list(bases), 'rational': rat is not None, 'reconstruction': rinfo,
              'seconds': round(time.time() - t0, 1), 'tag': a.tag}
    if rat is not None and a.numerators and rat:
        nums, ninfo = numerators_rational(F, A, rat, primes=core.PRIMES[:4], verbose=True)
        save = json.load(open(out_path))
        save['numerators_rational'] = None if nums is None else [[poly_to_json(c) for c in P] if P else None for P in nums]
        save['numerators_info'] = ninfo
        with open(out_path, 'w') as f: json.dump(save, f)
        logrec['numerators'] = nums is not None
    w2log(logrec)
    print(json.dumps(logrec), flush=True)

if __name__ == '__main__':
    main()
