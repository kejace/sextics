"""Compute E1-kernel reps at degree d, parity parD for several primes; CRT-reconstruct; save.
Usage: python3 w2_kernel.py d parD [nprimes]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

d = int(sys.argv[1]); parD = int(sys.argv[2]); npr = int(sys.argv[3]) if len(sys.argv) > 3 else 3
F, A, q, rec = load()
reps_by_prime, infos = [], []
for p in PRIMES[:npr]:
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    reps, info, red = kernel_reps(Fp, Ap, d, parD, p, verbose=True)
    reps_by_prime.append(reps); infos.append(info)
    # save mod-p reps
    with open(f'{RES}/w2/kernel_d{d}_p{parD}_mod{p}.json', 'w') as f:
        json.dump({'d': d, 'parD': parD, 'p': p, 'info': info, 'reps': [core.to_str(r) for r in reps]}, f)
    wlog(batch='e1-kernel', family='mult_kernel', degree=d, parity=parD, p=p, kernel_mod_J=info['kernel_mod_J'],
         dim_J=info['dim J_d'], dim_R=info['dim R_d'], rows=info['rows'], cols=info['cols'], seconds=info['seconds'])
    if len(set(i['kernel_mod_J'] for i in infos)) > 1:
        print('WARNING: kernel dimension differs between primes', infos, flush=True)
rat, bad = crt_reconstruct(reps_by_prime, PRIMES[:npr])
if rat is None:
    print(f'rational reconstruction failed (rep {bad}) with {npr} primes', flush=True)
    wlog(batch='e1-kernel-recon', degree=d, parity=parD, nprimes=npr, ok=False, failed_rep=bad)
else:
    # verify E1 for each over a fresh prime
    p3 = PRIMES[npr]
    out = []
    for k, D in enumerate(rat):
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p3)
        out.append({'d': d, 'parD': parD, 'D': core.to_str(D), 'nterms': len(D), 'E1_p_check': ok, 'maxden': max(abs(v.denominator) for v in D.values()), 'maxnum': max(abs(v.numerator) for v in D.values())})
        print(f'rep {k}: {len(D)} terms, E1 mod {p3}: {ok}, maxden {out[-1]["maxden"]}', flush=True)
    save_poly_json(f'{RES}/w2/kernel_d{d}_p{parD}_rational.json', out)
    wlog(batch='e1-kernel-recon', degree=d, parity=parD, nprimes=npr, ok=True, n=len(out), e1_checks=[o['E1_p_check'] for o in out])
print('done', flush=True)
