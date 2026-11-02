import sys, time
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
import fplib, numpy as np
F, A, q, G = fplib.load(); ch = fplib.Chars(G)
for d in [int(x) for x in sys.argv[1].split(',')]:
    irr = [c for c in sorted(ch.all_chars(d)) if fplib.is_irreducible_sector(ch, d, c)]
    print('d', d, 'irreducible sectors', irr, flush=True)
    for c in irr:
        t0 = time.time()
        B, Dm, info = fplib.sector_basis_numeric(F, A, q, d, c, ch)
        if info['adm'] == 0:
            print('  ', c, 'adm 0', info, flush=True); continue
        fam = fplib.SectorVarPro(F, A, q, B, Dm, c, ch, tag=f'd{d}')
        t1 = time.time()
        r = fam.minimize_from(1, maxiter=50)
        t2 = time.time()
        r2 = fam.lm_from(1, maxiter=10, init=r['u'])
        t3 = time.time()
        print('  ', c, 'adm', info['adm'], 'gap', info['sv_gap'], 'dimP', fam.dimP, 'rows', fam.R, 'nN', fam.nN, 'nQ', fam.nQ,
              'build', round(t1 - t0, 1), 'lbfgs50', round(t2 - t1, 1), 'rel', f"{r['relres']:.2e}", 'lm10', round(t3 - t2, 1), 'rel', f"{r2['relres']:.2e}", flush=True)
