"""Served Y-odd target space of a denominator D (XY frame) via compressed rows (feasible at d=8), then
sector analysis: classes in R_6 = R_f (x) R_f (U,V frame), bidegree-block projections, K3 block (6,0)+(0,6).
usage: leanimage.py '<label>' '<expr>' [uv|xy]"""
import sys, time, json, resource
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from flint import nmod_mat
from hodge import core
from hodge.fixed_d import FixedDSystem, _sparse_matrix
from hodge.core import monomials
from fastsolve import compress_chunked
from runD2 import uv_to_xy  # safe: runD2 guarded below

def served_space_lean(F, q, D, p, verbose=True):
    t0 = time.time()
    A0 = {monomials(6 * q, 1)[0]: 1}
    S = FixedDSystem(F, A0, q, D, p, with_target_space=True, verbose=verbose)
    M = _sparse_matrix(S, with_rhs=False)
    C = S.ncols; k = min(S.nrows, C + 40)
    Cm = compress_chunked(M, p, k); del M
    t1 = time.time()
    Mn = nmod_mat(k, C, Cm.ravel().tolist(), p); del Cm
    X, nullity = Mn.nullspace()
    t2 = time.time()
    a_idx = [j for j, c in enumerate(S.cols) if c[0] == 'A']
    if nullity == 0:
        return 0, [], dict(rows=S.nrows, cols=C, k=k, nullity=0, seconds=round(t2-t0,1))
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], p)
    Tr, rk = proj.transpose().rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[S.cols[j][2]] = v
        basis.append(poly)
    if verbose:
        print(f'[image d={S.d}] rows={S.nrows} cols={C} k={k} nullity={nullity} served_dim={rk} compress {t1-t0:.1f}s null {t2-t1:.1f}s', flush=True)
    return rk, basis, dict(rows=S.nrows, cols=C, k=k, nullity=nullity, served_dim=rk, seconds=round(t2-t0,1))

if __name__ == '__main__':
    label, expr = sys.argv[1], sys.argv[2]
    frame = sys.argv[3] if len(sys.argv) > 3 else 'uv'
    F, A, q, rec = core.load_instance('polarized-k3-remainder')
    D = core.parse(expr)
    if frame == 'uv': D = uv_to_xy(D)
    p = core.PRIMES[0]
    dim, basis, info = served_space_lean(F, q, D, p)
    out = {'label': label, 'expr': expr, 'frame': frame, 'D_xy': core.to_str(D), 'dim': dim, 'info': info,
           'basis': [core.to_str(b) for b in basis]}
    json.dump(out, open(f'target-image-{label}.json', 'w'))
    import k3proj
    an = k3proj.analyse(f'target-image-{label}.json')
    # is the instance target A in the served span (mod p)?
    T = k3proj.T
    vecs = [T.vector_of(k3proj.xy_to_uv_modp(core.parse(s, p)), 6) for s in out['basis']]
    vA = T.vector_of(k3proj.xy_to_uv_modp(core.reduce_mod(A, p)), 6)
    n = T.dim(6)
    if vecs:
        M1 = nmod_mat(len(vecs), n, [x for v in vecs for x in v], p).rank()
        M2 = nmod_mat(len(vecs)+1, n, [x for v in vecs for x in v] + vA, p).rank()
        an['target_in_span_mod_J'] = (M1 == M2)
    else:
        an['target_in_span_mod_J'] = False
    an['info'] = info; an['maxrss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    print(json.dumps(an), flush=True)
    with open('images.jsonl', 'a') as fh: fh.write(json.dumps(an) + '\n')
    from logutil import log
    log(label, core.degree(D), 'served-space', None, bool(an['target_in_span_mod_J']), -1, info.get('nullity', -1), info.get('seconds', 0),
        f"leanimage (consistency = target class in served span mod J); rows={info['rows']} cols={info['cols']} k={info['k']}; served_dim={dim}; class_dim={an.get('class_dim_in_R6', an.get('class_dim'))}; K3_proj={an.get('K3_projection_dim')}; blocks={an.get('block_projection_dims')}; frame={frame}; expr={expr}")
