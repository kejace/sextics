import sys, time
sys.path.insert(0, '.')
from hodge import core, fixed_d

for iid in ['collision-h3', 'hesse-triangle']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    t0 = time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    print(f'{iid}: e1_solvable={ok} {info} ({time.time()-t0:.1f}s)')
    t0 = time.time()
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True)
    print(f'{iid}: status={status} ', {k: v for k, v in data.items() if k not in ('N', 'Q')}, f'({time.time()-t0:.1f}s)')
    if status == 'solved':
        print('  N term counts', [len(n) for n in data['N']], 'Q terms', len(data['Q']))
    # negative control: perturb D
    D2 = dict(D); e = next(iter(D2)); D2[e] = D2[e] + 1
    t0 = time.time()
    status, data = fixed_d.solve_fixed_D(F, A, q, D2, max_primes=1)
    print(f'{iid}: perturbed D status={status}', data, f'({time.time()-t0:.1f}s)')
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F, q, D)
    print(f'{iid}: target image dim={dim} ({time.time()-t0:.1f}s)')
    for b in basis[:5]:
        print('   ', core.to_str({k: v if v < S.p//2 else v - S.p for k, v in b.items()})[:120])
