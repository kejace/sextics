"""Job: excl-parity<d>.  Parity-class closures at degree d: the 8 classes s = (b0, b2, y) of the (Z/2)^3 quotient of the
G_aff character group (all six h1-characters mixed inside a class, as in fqlx's degree-7 analysis).  Family = union of the
gauge-fixed E1-admissible families of the pure sectors in the class; Hessian relaxation quadrics; msolve -g 1 verdict
(fallback Singular); quads exported as JSON for branchkill.py.  Usage: parity_closure.py d [classes like 000,011] [threads]"""
import sys, time, json, os
from w2lib import *
from family import *
from msverdict import verdict
d = int(sys.argv[1])
want = sys.argv[2].split(',') if len(sys.argv) > 2 and sys.argv[2] != 'all' else None
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
OUT = PK + f'/parity_d{d}'
os.makedirs(OUT, exist_ok=True)
chars = sorted({chi_of({e: 1}) for e in monomials(d)})
classes = {}
for c in chars:
    classes.setdefault((c[0], c[1], c[3]), []).append(c)
for s, cs in sorted(classes.items()):
    tag = ''.join(map(str, s))
    if want and tag not in want: continue
    t0 = time.time()
    Ds, Ns, chis, gauge_idx = [], [], [], []
    per = {}
    for c in cs:
        D_, N_, ng, nj, ne = gauge_and_jacobian_family(d, c)
        base = len(Ds)
        gauge_idx.extend(range(base, base + ng))
        Ds.extend(D_); Ns.extend(N_); chis.extend([c] * len(D_))
        per[str(c)] = {'n': len(D_), 'gauge': ng, 'jac': nj, 'extra': ne}
    n = len(Ds)
    quads, qi = quadrics_from_family(Ds, Ns, d, verbose=False)
    keep = [k for k in range(n) if k not in set(gauge_idx)]
    kidx = {k: i for i, k in enumerate(keep)}
    quads2 = []
    for cf in quads:
        cf2 = {(kidx[k], kidx[l]): v for (k, l), v in cf.items() if k in kidx and l in kidx}
        if cf2: quads2.append(cf2)
    m = len(keep)
    rank, forced, pairs = quad_rank_and_squares(quads2, m)
    prefix = OUT + f'/parity_d{d}_{tag}'
    neq, nv = write_systems(quads2, m, prefix)
    json.dump({'p': p, 'chi': [list(chis[k]) for k in keep], 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads2],
               'D': [core.to_str(Ds[k]) for k in keep]}, open(prefix + '.quads.json', 'w'))
    info = {'d': d, 'parity_class': tag, 'sectors': per, 'n_unknowns': m, 'n_quadrics': len(quads2), 'quadric_rank': rank,
            'n_pair_monomials': len(pairs), 'forced_zero_by_squares': forced}
    if neq == 0:
        info['result'] = 'NO EQUATIONS'
    else:
        v = verdict(prefix + '.ms', threads=threads, timeout=int(os.environ.get('MSOLVE_TIMEOUT', '7200')), verbose=1, logfile=prefix + '.msolve.log')
        info['msolve'] = {k: val for k, val in v.items() if k not in ('log_tail', 'pure_powers')}
        info['result'] = 'EXCLUDED (only c=0)' if v.get('verdict', '').startswith('ONLY-ORIGIN') else v.get('verdict')
        if not v.get('verdict', '').startswith('ONLY-ORIGIN') and m <= 60:
            sv = singular_verdict(prefix + '.sing', timeout=int(os.environ.get('SING_TIMEOUT', '7200')))
            info['singular'] = {k: sv[k] for k in ('sdim', 'svdim', 'rc', 'seconds')}
            if sv['sdim'] == 0: info['result'] = 'EXCLUDED (only c=0; Singular)'
    info['seconds'] = round(time.time() - t0, 1)
    log('closure-parity-class', **info, note='Hessian relaxation of the full E1-admissible family with all h1-characters mixed inside the parity class (b0,b2,y)')
print('DONE', file=sys.stderr)
