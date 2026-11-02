"""(c)/(b): sector-wise exhaustive Hessian relaxation for SEMI-INVARIANT denominators.
For each degree d and each G-character sector: E1 space V (incl. J_d part), quadrics of the
relaxation E2 mod J in dim V variables, Groebner basis in Singular (char p), dimension/degree.
dim 0 (affine) => projective variety empty => no semi-invariant D in this sector (rigorous via the
Hilbert-polynomial comparison, see commentary).  dim 1 => finitely many candidate lines (listed).
"""
import sys, json, time, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
p = P2
degs = [int(x) for x in sys.argv[1:]] or [6, 7, 8]
JS = JacobianSectors(F2, p)
kA = key_of_poly(A, JS.gens)
OUT = '/tmp/hodge/results/cycle-q2/wave2-excl/receipts'
for d in degs:
    t0 = time.time()
    summary = {}
    degB = 2 * d + 6
    for key in all_sector_keys(JS.gens, d):
        r = JS.e1_kernel(A, d, key)
        k = r['kernel_dim']
        if k == 0:
            continue
        Ds = r['basis']
        N0, B = relaxation_data(JS, A, Ds)
        quadrics = {}
        for (a, b), Bab in B.items():
            if not Bab:
                continue
            psi = key_of_poly(Bab, JS.gens)
            C = JS.normal_form_coords([Bab], degB, psi)
            for s in range(C.ncols()):
                v = int(C[0, s])
                if v:
                    quadrics.setdefault((psi, s), {})[(a, b)] = v
        vars_ = [f'c{i}' for i in range(k)]
        polys = []
        for (psi, s), terms in quadrics.items():
            mon = []
            for (a, b), v in terms.items():
                coef = v if a == b else (2 * v) % p
                mon.append(f'{coef}*c{a}*c{b}' if a != b else f'{coef}*c{a}^2')
            polys.append('+'.join(mon))
        tag = f'd{d}_k' + '_'.join(map(str, key))
        sing = f'ring r = {p},({",".join(vars_)}),dp;\n'
        sing += 'ideal I = ' + (',\n'.join(polys) if polys else '0') + ';\n'
        sing += 'option(redSB); ideal G = std(I);\nprint("DIM " + string(dim(G)));\nprint("MULT " + string(mult(G)));\nprint("SIZE " + string(size(G)));\nprint("GB"); print(G); quit;\n'
        fn = f'{OUT}/relax_{tag}.sing'
        open(fn, 'w').write(sing)
        out = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=3600).stdout
        dim_ = int([l for l in out.splitlines() if l.startswith('DIM')][0].split()[1])
        mult_ = int([l for l in out.splitlines() if l.startswith('MULT')][0].split()[1])
        gb = out.split('GB\n', 1)[1].strip() if 'GB\n' in out else ''
        summary[','.join(map(str, key))] = dict(kernel_dim=k, ann_dim=r['dim'], jdim=r['jdim'], n_quadrics=len(quadrics), dim=dim_, mult=mult_,
                                                basis=[to_str(Dk) for Dk in Ds], gb=gb[:4000])
        print(d, key, 'k', k, 'ann', r['dim'], 'J', r['jdim'], 'quadrics', len(quadrics), 'DIM', dim_, 'MULT', mult_, flush=True)
    json.dump(summary, open(f'{OUT}/sector_relax_d{d}.json', 'w'), indent=1)
    nsec = len(summary); empty = sum(1 for v in summary.values() if v['dim'] == 0); pts = sum(1 for v in summary.values() if v['dim'] == 1)
    pos = sum(1 for v in summary.values() if v['dim'] >= 2)
    rec = dict(exp='C-sector-relaxation', engine=ENGINE_TAG, d=d, p=p, sectors_with_E1=nsec, sectors_empty=empty, sectors_finite=pts, sectors_posdim=pos,
               total_e1_dim=sum(v['kernel_dim'] for v in summary.values()), seconds=round(time.time() - t0, 1), file=f'{OUT}/sector_relax_d{d}.json',
               finite=[(k, v['mult'], v['kernel_dim']) for k, v in summary.items() if v['dim'] == 1],
               posdim=[(k, v['dim'], v['mult'], v['kernel_dim']) for k, v in summary.items() if v['dim'] >= 2])
    print(json.dumps(rec)); xlib.log(rec)
