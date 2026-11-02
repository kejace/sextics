"""Structure of the degree-7 E1 kernel: restrict E1 to the base-type space B_7 = span{P(X) M^a g^b} (and the mirror
space span{Q(Y) L^a f^b}), find kernel mod J, locate it in the (a,b)-blocks and reconstruct rational generators."""
import sys, json, time, itertools
import lib
from lib import *
from flint import nmod_mat
p1, p2 = core.PRIMES[0], core.PRIMES[1]
TARGET = sys.argv[1] if len(sys.argv) > 1 else 'X0**3*Y0**3'
Tq = core.parse(TARGET); degT = core.degree(Tq); parT_ = core.parity_of(Tq)
d = int(sys.argv[2]) if len(sys.argv) > 2 else 7
SIDE = sys.argv[3] if len(sys.argv) > 3 else 'Y'   # 'Y': P(X) M^a g^b ; 'X': Q(Y) L^a f^b

def xmons(k, ys=False):
    if ys: return [m for m in core.monomials(k) if m[0] == m[1] == m[2] == 0]
    return [m for m in core.monomials(k) if m[3] == m[4] == m[5] == 0]

def blocks(d):
    """list of (label, [polys]) of base type and total degree d."""
    out = []
    Q2, Q6 = (M, g) if SIDE == 'Y' else (L, f)
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0: continue
            W = core.power(Q2, a) if a else {(0,)*6: Fraction(1)}
            if b: W = core.mul(W, Q6)
            polys = [core.mul({m: Fraction(1)}, W) for m in xmons(k, ys=(SIDE == 'X'))]
            if polys: out.append((f'S_{k}*{"M" if SIDE=="Y" else "L"}^{a}*{"g" if SIDE=="Y" else "f"}^{b}', polys))
    return out

def rowspace_mod(polys, mons, p):
    idx = {m: k for k, m in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for j, c in enumerate(polys):
        for e, v in core.reduce_mod(c, p).items(): flat[j * len(mons) + idx[e]] = v
    return nmod_mat(len(polys), len(mons), flat, p)

def kernel_coords(T, Bpolys, p):
    """coords (in Bpolys basis) of {D in span(B): T*D in J} and of {D in span(B): D in J}; returns (E1rows, Jrows) as nmod_mat rref."""
    Fp = core.reduce_mod(F, p); Tp = core.reduce_mod(T, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    nB = len(Bpolys)
    def space(mult, deg_tot):
        # columns: J_{deg_tot} generators (all parities) and mult*B_j ; nullspace projected onto B-coords
        rows = core.monomials(deg_tot); idx = {m: k for k, m in enumerate(rows)}
        cols = []
        for i in range(6):
            for m in core.monomials(deg_tot - 5):
                cols.append(core.shift(Fi[i], m, 1, p))
        nJ = len(cols)
        for b in Bpolys:
            cols.append(core.mul(core.reduce_mod(b, p), mult, p) if mult is not None else core.reduce_mod(b, p))
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, c in enumerate(cols):
            for e, v in c.items(): flat[idx[e] * C + j] = v
        X, nul = nmod_mat(R, C, flat, p).nullspace()
        if nul == 0: return nmod_mat(0, nB, [], p), 0
        proj = nmod_mat(nB, nul, [int(X[nJ + j, k]) for j in range(nB) for k in range(nul)], p)
        Tr, rk = proj.transpose().rref()
        return nmod_mat(rk, nB, [int(Tr[i, j]) for i in range(rk) for j in range(nB)], p), rk
    E1, rkE = space(Tp, d + degT)
    Jr, rkJ = space(None, d)
    return E1, rkE, Jr, rkJ

def pivots(Mr, rk, n):
    out = []
    for i in range(rk):
        j = 0
        while j < n and int(Mr[i, j]) == 0: j += 1
        out.append(j)
    return out

def canon_complement(E1, rkE, Jr, rkJ, nB, p):
    """rows of rref(E1) whose pivots are not J-pivots, after full reduction against J rows -> canonical kernel reps mod J."""
    rowsAll = [[int(Jr[i, j]) for j in range(nB)] for i in range(rkJ)] + [[int(E1[i, j]) for j in range(nB)] for i in range(rkE)]
    Sr, rkS = nmod_mat(len(rowsAll), nB, [v for r in rowsAll for v in r], p).rref()
    assert rkS == rkE, (rkS, rkE)
    pivJ = set(pivots(Jr, rkJ, nB))
    ker = []
    for i in range(rkS):
        j = 0
        while j < nB and int(Sr[i, j]) == 0: j += 1
        if j not in pivJ: ker.append([int(Sr[i, jj]) for jj in range(nB)])
    return ker

if __name__ == '__main__':
    t0 = time.time()
    bl = blocks(d)
    labels = []; Bpolys = []
    for lab, polys in bl:
        for pp_ in polys:
            labels.append(lab); Bpolys.append(pp_)
    nB = len(Bpolys)
    print(f'TARGET={TARGET} d={d} SIDE={SIDE}: base-type space dim {nB}; blocks: {[(lab, len(pl)) for lab, pl in bl]}', flush=True)
    res = {}
    for p in (p1, p2):
        E1, rkE, Jr, rkJ = kernel_coords(Tq, Bpolys, p)
        ker = canon_complement(E1, rkE, Jr, rkJ, nB, p)
        print(f'  p={p}: dim(E1 cap B)={rkE}, dim(J cap B)={rkJ}, kernel mod J = {len(ker)}  ({time.time()-t0:.1f}s)', flush=True)
        res[p] = (E1, rkE, Jr, rkJ, ker)
        # block localisation: for each union of blocks (prefix-closed subsets), kernel dim within
        if p == p1:
            blab = [lab for lab, _ in bl]
            for r in range(1, len(blab) + 1):
                for sub in itertools.combinations(range(len(blab)), r):
                    cols_keep = [j for j, lab in enumerate(labels) if lab in {blab[s] for s in sub}]
                    # E1 cap span(sub): rows of E1 restricted... need intersection: solve E1-rowspace cap coordinate subspace
                    def inter_dim(Mr, rk):
                        if rk == 0: return 0
                        # vectors in rowspace with zero outside cols_keep: nullspace of (coeff^T restricted to complement cols)
                        comp = [j for j in range(nB) if j not in set(cols_keep)]
                        if not comp: return rk
                        sub_m = nmod_mat(len(comp), rk, [int(Mr[i, j]) for j in comp for i in range(rk)], p)
                        _, nul = sub_m.nullspace()
                        return nul
                    dE = inter_dim(E1, rkE); dJ = inter_dim(Jr, rkJ)
                    if dE - dJ > 0:
                        print(f'    blocks {[blab[s] for s in sub]}: dim E1-cap={dE}, J-cap={dJ}, kernel mod J = {dE - dJ}')
    # rational reconstruction of canonical kernel reps (same pivot structure assumed)
    k1 = res[p1][4]; k2 = res[p2][4]
    gens = []
    for v1, v2 in zip(k1, k2):
        poly = {}; coords = []
        ok = True
        for j in range(nB):
            a, mm = core.crt_pair(v1[j], p1, v2[j], p2)
            r = core.ratrecon(a, mm)
            if r is None: ok = False; break
            if r: coords.append((labels[j], j, r)); poly = core.add(poly, core.scale(Bpolys[j], r))
        if ok:
            print('  generator (structured coords):', [(lab, str(r)) for lab, j, r in coords][:40], '... nterms', len(coords))
            gens.append(core.to_str(poly))
        else:
            print('  ratrecon failed for a generator')
    json.dump({'target': TARGET, 'd': d, 'side': SIDE, 'labels': labels, 'generators': gens}, open(f'{WORK}/kstruct_{TARGET.replace("*","")}_d{d}_{SIDE}.json', 'w'))
    print('done', time.time() - t0)
