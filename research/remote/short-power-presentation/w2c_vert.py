"""Job: restricted E1 kernels of [A] on (a) C[X0,X1,X2,M]_d, (b) C[L,Y0,Y1,Y2]_d, (c) their sum U_d (an independent
basis is extracted first; U_d contains the vertical ring C[e(X),e(Y),L,M]_d for d<=11), (d) the vertical ring itself.
Usage: python3 w2c_vert.py <out.jsonl> d1 d2 ...   (both parities each)."""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2c_lib import *
import flint
try: flint.ctx.threads = int(os.environ.get('FLINT_THREADS', '8'))
except Exception: pass

def independent_subfamily(fam, p):
    """Select a linearly independent subset of the family polynomials (mod p)."""
    polys = [reduce_mod(P, p) for _, P in fam]
    rows = {}
    for P in polys:
        for e in P: rows.setdefault(e, len(rows))
    R, C = len(rows), len(polys)
    flat = [0] * (R * C)
    for j, P in enumerate(polys):
        for e, v in P.items(): flat[rows[e] * C + j] = v
    Mx = nmod_mat(R, C, flat, p)
    Mr, rk = Mx.rref()
    piv = []
    for i in range(rk):
        j = 0
        while j < C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    return [fam[j] for j in piv]

out_path = sys.argv[1]
degs = [int(x) for x in sys.argv[2:]]
F, A, q, rec = core.load_instance(IID)
p = core.PRIMES[0]
for d in degs:
    for par in (0, 1):
        famXM = fam_coneXM(d, par); famLY = fam_coneLY(d, par); famV = fam_vert(d, par)
        famU = independent_subfamily(famXM + famLY, p)
        res = {}
        t0 = time.time()
        for fname, fam in (('coneXM', famXM), ('coneLY', famLY), ('union', famU), ('vert', famV)):
            if not fam: res[fname] = {'fam_dim': 0, 'kernel_dim': 0, 'kernel': []}; continue
            vecs, info = restricted_kernel(F, A, fam, d, p)
            an = analyze_kernel(F, vecs, fam, d, p)
            res[fname] = {'fam_dim': info['fam_dim'], 'rows': info['rows'], 'cols': info['cols'], 'kernel_dim': info['kernel_dim'],
                          'F_multiples': sum(1 for a in an if a['F_divides']), 'in_J': sum(1 for a in an if a['in_J']),
                          'kernel': an if info['kernel_dim'] <= 12 else 'omitted (%d)' % info['kernel_dim']}
        rec_out = {'exp': 'restricted-E1-kernel', 'd': d, 'parD': par, 'prime': p, 'families': res,
                   'seconds': round(time.time() - t0, 1), 'job': 'construct-spp-vert',
                   'result': ' | '.join('%s: fam %d kernel %d (Fmult %d, inJ %d)' % (k, v['fam_dim'], v['kernel_dim'], v.get('F_multiples', 0), v.get('in_J', 0)) for k, v in res.items())}
        with open(out_path, 'a') as f: f.write(json.dumps(rec_out) + '\n')
        print('DONE', d, par, rec_out['result'], flush=True)
print('ALL DONE', flush=True)
