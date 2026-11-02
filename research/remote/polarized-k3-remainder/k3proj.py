"""Sector decomposition of a served-target space: convert target_image basis (XY frame, mod P) to the (U,V)
frame, take classes in R_6 = (+)_{a+b=6} R_a (x) R_b, and report the dimension of the projection onto each
bidegree block, in particular the K3 block (6,0)+(0,6)."""
import sys, json, glob; sys.path.insert(0, '.')
from uvring import *
from flint import nmod_mat
f = load_f(); R = TernaryJacobianRing(f, maxdeg=12); T = TensorRing(R, maxdeg=12)

def xy_to_uv_modp(Pxy, p=P):
    """X = (U+V)/2, Y = (U-V)/2; slots: X->U(0..2), Y->V(3..5)."""
    inv2 = pow(2, -1, p)
    def var(i): e=[0]*6; e[i]=1; return {tuple(e): 1}
    U = [var(i) for i in range(3)]; V = [var(3+i) for i in range(3)]
    imgs = [core.scale(core.add(U[i], V[i], p), inv2, p) for i in range(3)] + [core.scale(core.sub(U[i], V[i], p), inv2, p) for i in range(3)]
    out = {}
    for e, c in Pxy.items():
        term = {(0,)*6: c % p}
        for i in range(6):
            if e[i]: term = core.mul(term, core.power(imgs[i], e[i], p), p)
        out = core.add(out, term, p)
    return out

def analyse(path):
    J = json.load(open(path))
    basis = [core.parse(s, P) for s in J['basis']]
    vecs = [T.vector_of(xy_to_uv_modp(b), 6) for b in basis]
    n = T.dim(6)
    out = {'label': J['label'], 'served_dim_XY': J['dim']}
    if not vecs:
        out['class_dim'] = 0; return out
    M = nmod_mat(len(vecs), n, [x for v in vecs for x in v], P)
    out['class_dim_in_R6'] = M.rank()
    blocks = {}
    for (a, b) in [(6,0),(5,1),(4,2),(3,3),(2,4),(1,5),(0,6)]:
        cols = [i for i, t in enumerate(T.bases[6]) if (t[0], t[1]) == (a, b)]
        if not cols: blocks[f'({a},{b})'] = 0; continue
        Mb = nmod_mat(len(vecs), len(cols), [v[i] for v in vecs for i in cols], P)
        blocks[f'({a},{b})'] = Mb.rank()
    out['block_projection_dims'] = blocks
    colsK3 = [i for i, t in enumerate(T.bases[6]) if t[0] in (0, 6)]
    MK = nmod_mat(len(vecs), len(colsK3), [v[i] for v in vecs for i in colsK3], P)
    out['K3_projection_dim'] = MK.rank()
    # residue character sectors by V-degree b mod 6: b in {0,6} K3; others supported
    return out

if __name__ == '__main__':
    files = sys.argv[1:] or sorted(glob.glob('target-image-*.json'))
    res = []
    for fpath in files:
        r = analyse(fpath); res.append(r); print(r, flush=True)
    json.dump(res, open('k3-projections.json', 'w'), indent=1)
