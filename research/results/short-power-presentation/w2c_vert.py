"""Job: restricted E1 kernels of [A] on the union family U_d = C[X0,X1,X2,M]_d + C[L,Y0,Y1,Y2]_d (sums of a cone with
vertex the line {X=M=0} and a cone with vertex {Y=L=0}; contains the 'vertical' ring C[e(X),e(Y),L,M]_d for d<=11),
one nullspace per (d, parity); then the dimensions of the kernel's intersections with the sub-families.
Usage: python3 w2c_vert.py <out.jsonl> d1 d2 ...   (both parities each)."""
import sys, json, time, os
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2c_lib import *
import flint
try: flint.ctx.threads = int(os.environ.get('FLINT_THREADS', '8'))
except Exception: pass

def sub_dim(vecs, names_sub, p):
    """dimension of span(vecs) ∩ span(sub-family coordinates): vectors supported only on names_sub."""
    if not vecs: return 0
    names = sorted({n for v in vecs for n in v})
    idx = {n: i for i, n in enumerate(names)}
    other = [i for i, n in enumerate(names) if n not in names_sub]
    K = len(vecs)
    # kernel vectors c with sum_k c_k v_k supported on names_sub: solve the 'other' coordinates = 0
    rows = [[v.get(names[i], 0) for v in vecs] for i in other]
    if not rows: return K
    Mx = nmod_mat(len(rows), K, [x for r in rows for x in r], p)
    return K - Mx.rank()

out_path = sys.argv[1]
degs = [int(x) for x in sys.argv[2:]]
F, A, q, rec = core.load_instance(IID)
p = core.PRIMES[0]
for d in degs:
    for par in (0, 1):
        famXM = fam_coneXM(d, par); famLY = fam_coneLY(d, par); famV = fam_vert(d, par)
        # union: XM ∪ LY, dropping duplicates of pure C[L,M] monomials (present in both)
        seen = set(); fam = []
        for n, P in famXM + famLY:
            key = tuple(sorted(P.items()))
            if key in seen: continue
            seen.add(key); fam.append((n, P))
        t0 = time.time()
        vecs, info = restricted_kernel(F, A, fam, d, p)
        an = analyze_kernel(F, vecs, fam, d, p)
        namesXM = {n for n, _ in famXM}; namesLY = {n for n, _ in famLY}
        dXM = sub_dim(vecs, namesXM, p); dLY = sub_dim(vecs, namesLY, p)
        # vertical sub-family: express fam_vert elements in the union basis is nontrivial; instead test membership of each vert element's kernel separately
        vV, infoV = restricted_kernel(F, A, famV, d, p) if famV else ([], {'kernel_dim': 0})
        anV = analyze_kernel(F, vV, famV, d, p)
        rec_out = {'exp': 'restricted-E1-kernel', 'family': 'XM+YL (union of cones)', 'd': d, 'parD': par, 'prime': p,
                   'fam_dim': info['fam_dim'], 'rows': info['rows'], 'cols': info['cols'], 'kernel_dim': info['kernel_dim'],
                   'kernel_in_coneXM': dXM, 'kernel_in_coneLY': dLY, 'kernel_F_multiples': sum(1 for a in an if a['F_divides']),
                   'kernel_in_J': sum(1 for a in an if a['in_J']),
                   'vert_fam_dim': len(famV), 'vert_kernel_dim': infoV['kernel_dim'], 'vert_kernel': anV,
                   'seconds': round(time.time() - t0, 1), 'job': 'construct-spp-vert',
                   'result': 'union kernel=%d (XM %d, LY %d, F-mult %d, inJ %d); vert kernel=%d' % (info['kernel_dim'], dXM, dLY, sum(1 for a in an if a['F_divides']), sum(1 for a in an if a['in_J']), infoV['kernel_dim'])}
        if info['kernel_dim'] <= 12:
            rec_out['kernel'] = an
        with open(out_path, 'a') as f: f.write(json.dumps(rec_out) + '\n')
        print('DONE', d, par, rec_out['result'], flush=True)
print('ALL DONE', flush=True)
