"""Step 1: E1-admissible denominator spaces per degree d=4..12 and Y-parity (the only character of F):
dim (R/J)_d, dim J_d, non-Jacobian kernel of [A] mod J, total E1-admissible dim, gauge F*R_{d-6}."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
import w2k3 as W
from hodge import core, jacobian
F, A, q = W.load()
d0, d1 = int(sys.argv[1]), int(sys.argv[2])
pars = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0, 1]
for d in range(d0, d1 + 1):
    for par in pars:
        t0 = time.time()
        mk = jacobian.mult_kernel(F, A, d, parity=par)
        gauge = len(core.monomials(d - 6, par)) if d >= 6 else 0
        rec = {'step': 'dims', 'd': d, 'parity': par, 'sector': 'Y-even' if par == 0 else 'Y-odd',
               'ring_dim': mk['ring_dim'], 'jacobian_dim': mk['jacobian_dim'], 'nonjac_kernel': mk['dim'],
               'E1_admissible_dim': mk['kernel_dim'] if mk['kernel_dim'] else mk['jacobian_dim'] + mk['dim'],
               'gauge': gauge, 'p': mk['p'], 'rows': mk['rows'], 'cols': mk['cols'], 'seconds': round(time.time() - t0, 1),
               'note': 'mult_kernel: kernel of [A]:(R/J)_d -> (R/J)_{d+6}; E1-admissible = J_d + kernel reps; gauge = F*R_{d-6} (illegal directions)'}
        # mult_kernel's kernel_dim is dim{D in R_d: AD in J} before quotient (0 when computed from projection only)
        rec['E1_admissible_dim'] = mk['jacobian_dim'] + mk['dim']
        W.log_batch(rec)
        print(json.dumps(rec), flush=True)
        if d >= 10 and par == 0 and mk['dim']:
            json.dump({'d': d, 'basis': [{str(k): v for k, v in b.items()} for b in mk['basis']], 'p': mk['p']},
                      open(f'wave2-e1kernel-d{d}-even.json', 'w'))
