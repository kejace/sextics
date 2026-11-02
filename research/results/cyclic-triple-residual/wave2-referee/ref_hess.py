"""R3: Hessian relaxation  A * G^T Hess(F) G in J_{2d}  for D = G(F), own quadrics, fresh prime, own gauge slice (X0).
d = 5: exact over Q (Singular char 0: dim(std(I)) == 0) and mod p;  d = 6: mod p, both slices, Singular dim + msolve pure powers;
d = 7: all 18 character sectors (G-semi-invariant D), msolve pure powers (+ Singular dim)."""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-referee')
import ref_lib as L
P = int(os.environ.get('REF_P', '999999937'))
OUT = '/tmp/hodge/results/cyclic-triple-residual/wave2-referee'
LOG = OUT + '/referee.log.jsonl'
def log(rec):
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['script'] = 'ref_hess.py'
    with open(LOG, 'a') as f: f.write(json.dumps(rec, default=str) + '\n')
L.check_setup()

def test(d, p, slice_var, sector=None, tag='', singular=True, msolve=True, threads=2, timeout=1800, parity=None):
    t0 = time.time()
    th, quads, info = L.relaxation_quadrics(d, p, sector=sector, verbose=True, parity=parity)
    removed = L.gauge_slice(th, slice_var)
    keep, gq = L.restrict(quads, len(th), removed)
    names = [f'u{i}' for i in range(len(keep))]
    base = f'{OUT}/hess_d{d}_{"Q" if p is None else "p"+str(p)}_{slice_var}{tag}'
    rec = {'check': f'hessian_relaxation_d{d}{tag}', 'd': d, 'sector': str(sector), 'parity': parity, 'mode': 'exact Q' if p is None else f'mod {p}', 'slice': slice_var,
           'n_theta': len(th), 'n_removed': len(removed), 'n_vars': len(keep), 'n_quadrics': len(gq), 'rankJ_2d': info['rankJ_2d'], 'ci_dimJ_2d': info['ci_dimJ_2d'],
           'ci_ok': info['rankJ_2d'] == info['ci_dimJ_2d']}
    if not keep:
        rec['m_primary'] = True; rec['note'] = 'no unknowns after gauge fixing'
    else:
        polys = [L.quad_str(q, names, p) for q in gq]
        if singular:
            L.write_singular(base + '.sing', names, 0 if p is None else p, polys)
            rc, dt, out, err = L.run_singular(base + '.sing', timeout=timeout)
            rec['singular'] = {'rc': rc, 'dim': out[0] if out else None, 'vdim': out[1] if len(out) > 1 else None, 'seconds': round(dt, 1)}
            rec['m_primary_singular'] = (rc == 0 and bool(out) and out[0].strip() == '0')
        if msolve and p is not None:
            L.write_msolve(base + '.ms', names, p, polys)
            rc, dt, lead = L.msolve_lead(base + '.ms', base + '_lead.txt', threads=threads, timeout=timeout)
            pure, allp = L.pure_powers(lead, names)
            rec['msolve'] = {'rc': rc, 'lead_size': len(lead), 'pure_powers': len(pure), 'seconds': round(dt, 1)}
            rec['m_primary_msolve'] = (rc == 0 and allp)
    rec['seconds'] = round(time.time() - t0, 1)
    print('  ->', {k: v for k, v in rec.items() if k not in ('check',)}, flush=True)
    log(rec)
    return rec

which = sys.argv[1] if len(sys.argv) > 1 else 'all'
if which in ('all', 'd5'):
    test(5, None, 'X0', singular=True, msolve=False)          # exact over Q, Singular char 0
    test(5, P, 'X0')
    test(5, P, 'Y2')
if which in ('all', 'd6'):
    test(6, P, 'X0')
    test(6, P, 'Y2')
if which in ('all', 'd7'):
    for c in L.SECTORS:
        test(7, P, 'X0', sector=c, tag=f'_s{c[0]}{c[1]}', singular=True, msolve=True)
if which in ('d7yeven', 'd7yodd'):
    test(7, P, 'X0', parity=0 if which == 'd7yeven' else 1, tag='_' + which[2:], singular=False, msolve=True,
         threads=int(os.environ.get('REF_T', '4')), timeout=int(os.environ.get('REF_TO', '3000')))
if which == 'd7full':
    test(7, P, 'X0', singular=False, msolve=True, threads=int(os.environ.get('REF_T', '6')), timeout=int(os.environ.get('REF_TO', '1500')))
print('done', flush=True)
