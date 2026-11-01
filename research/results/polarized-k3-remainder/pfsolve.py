"""Parity-free fixed-D solver / served-target space in the (U,V) frame (slots 0..2 = U, 3..5 = V), mod P.
Builds the linear system  -q D N(F) + F (D div N - N(D)) - F^(q+1) Q + A' D^2 = 0  with ALL monomials
(no parity restriction) and target columns A' running over ALL degree-6q monomials, and returns the served
target space (as a list of polys) plus its dimension modulo J(F) via the tensor Jacobian ring."""
import sys, time, json; sys.path.insert(0, '.')
from uvring import *
from hodge.core import monomials, shift, mul, add, sub, scale, pderiv, power, NV
from flint import nmod_mat
p = P
uvf = json.load(open('uv_forms.json'))
Fh = core.parse(uvf['Fhat'], p); Ah = core.parse(uvf['Ahat'], p)
q = 1

def served_space(D, verbose=True):
    t0 = time.time()
    d = core.degree(D)
    degN = d + 6*q - 5; degQ = 2*d - 6
    Fi = [pderiv(Fh, i, p) for i in range(NV)]
    Di = [pderiv(D, i, p) for i in range(NV)]
    FD = mul(Fh, D, p)
    G = [scale(add(scale(mul(D, Fi[i], p), q, p), mul(Fh, Di[i], p), p), p-1, p) for i in range(NV)]
    Fq1 = scale(power(Fh, q+1, p), p-1, p)
    D2 = mul(D, D, p)
    cols, colpolys = [], []
    for i in range(NV):
        for m in monomials(degN):
            col = shift(G[i], m, 1, p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = add(col, shift(FD, tuple(m2), m[i], p), p)
            cols.append(('N', i, m)); colpolys.append(col)
    if degQ >= 0:
        for m in monomials(degQ):
            cols.append(('Q', -1, m)); colpolys.append(shift(Fq1, m, 1, p))
    tmons = monomials(6*q)
    for m in tmons:
        cols.append(('A', -1, m)); colpolys.append(shift(D2, m, 1, p))
    rows = {}
    for col in colpolys:
        for e in col:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    flat = [0] * (R * C)
    for j, col in enumerate(colpolys):
        for e, v in col.items():
            flat[rows[e] * C + j] = v
    M = nmod_mat(R, C, flat, p)
    t1 = time.time()
    X, nullity = M.nullspace()
    t2 = time.time()
    a_idx = [j for j, c in enumerate(cols) if c[0] == 'A']
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, k]) for j in a_idx for k in range(nullity)], p)
    T_ = proj.transpose(); Tr, rk = T_.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, j in enumerate(a_idx):
            v = int(Tr[i, jj])
            if v: poly[cols[j][2]] = v
        basis.append(poly)
    if verbose:
        print(f'[pf d={d}] rows={R} cols={C} (targets {len(a_idx)}) nullity={nullity} served_dim={rk} build {t1-t0:.1f}s null {t2-t1:.1f}s', flush=True)
    return basis, dict(d=d, rows=R, cols=C, nullity=nullity, served_dim=rk, seconds=round(t2-t0,1))

_ring = None
def ring():
    global _ring
    if _ring is None:
        f = load_f(); Rr = TernaryJacobianRing(f, maxdeg=12); _ring = TensorRing(Rr, maxdeg=12)
    return _ring

def class_analysis(basis):
    T = ring(); n = T.dim(6)
    vecs = [T.vector_of(b, 6) for b in basis]
    out = {}
    if not vecs:
        return dict(class_dim=0)
    M = nmod_mat(len(vecs), n, [x for v in vecs for x in v], p)
    out['class_dim_mod_J'] = M.rank()
    # swap-antisymmetric (tau-invariant residue) part: project classes v -> v - swap(v)
    Sw = swap_matrix(T, 6)
    Vm = nmod_mat(n, len(vecs), [v[i] for i in range(n) for v in vecs], p)
    SV = Sw * Vm
    anti = nmod_mat(n, len(vecs), [(int(Vm[i, j]) - int(SV[i, j])) % p for i in range(n) for j in range(len(vecs))], p)
    symm = nmod_mat(n, len(vecs), [(int(Vm[i, j]) + int(SV[i, j])) % p for i in range(n) for j in range(len(vecs))], p)
    out['class_dim_tau_invariant(anti-swap)'] = anti.rank(); out['class_dim_tau_anti(sym-swap)'] = symm.rank()
    blocks = {}
    for (a, b) in [(6,0),(5,1),(4,2),(3,3),(2,4),(1,5),(0,6)]:
        colsb = [i for i, t in enumerate(T.bases[6]) if (t[0], t[1]) == (a, b)]
        Mb = nmod_mat(len(vecs), len(colsb), [v[i] for v in vecs for i in colsb], p)
        blocks[f'({a},{b})'] = Mb.rank()
    out['block_projection_dims'] = blocks
    colsK3 = [i for i, t in enumerate(T.bases[6]) if t[0] in (0, 6)]
    MK = nmod_mat(len(vecs), len(colsK3), [v[i] for v in vecs for i in colsK3], p)
    out['K3_projection_dim'] = MK.rank()
    # is the instance target served?  class of Ahat in span?
    vA = T.vector_of(Ah, 6)
    M2 = nmod_mat(len(vecs) + 1, n, [x for v in vecs for x in v] + vA, p)
    out['target_in_span_mod_J'] = (M2.rank() == M.rank())
    return out

def sum_of_spaces(b1, b2):
    T = ring(); n = T.dim(6)
    vecs = [T.vector_of(b, 6) for b in b1 + b2]
    if not vecs: return 0
    return nmod_mat(len(vecs), n, [x for v in vecs for x in v], p).rank()

if __name__ == '__main__':
    label, expr = sys.argv[1], sys.argv[2]
    D = core.parse(expr, p)
    basis, info = served_space(D)
    an = class_analysis(basis)
    print(label, info, an, flush=True)
    json.dump({'label': label, 'expr': expr, 'info': info, 'analysis': an, 'basis': [core.to_str({e: v for e, v in b.items()}) for b in basis]},
              open(f'pf-{label}.json', 'w'))
