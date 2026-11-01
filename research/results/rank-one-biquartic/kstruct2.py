"""Parity-aware version: kernel of E1 restricted to the Y-shape base-type space B_d = span{P(X) M^a g^b}, modulo J_d
(and optionally modulo lower-degree kernel times linear forms), block localisation, rational generators via 2 primes."""
import sys, json, time, itertools
import lib
from lib import *
from flint import nmod_mat
p1, p2 = core.PRIMES[0], core.PRIMES[1]
TARGET = sys.argv[1] if len(sys.argv) > 1 else 'X0**3*Y0**3'
Tq = core.parse(TARGET); degT = core.degree(Tq); parT_ = core.parity_of(Tq)
d = int(sys.argv[2]) if len(sys.argv) > 2 else 7
parD = 0

def xmons(k):
    return [m for m in core.monomials(k) if m[3] == m[4] == m[5] == 0]

def blocks(d):
    out = []
    for a in range(0, d // 2 + 1):
        for b in range(0, 2):
            k = d - 2 * a - 6 * b
            if k < 0: continue
            W = core.power(M, a) if a else {(0,)*6: Fraction(1)}
            if b: W = core.mul(W, g)
            polys = [core.mul({m: Fraction(1)}, W) for m in xmons(k)]
            if polys: out.append((f'S{k}.M{a}.g{b}', polys))
    return out

def kernel_coords(T, Bpolys, p, extra=None):
    """rref row-matrices (coords in Bpolys basis) of E1 = {D in span B : T*D in J} and J' = {D in span B : D in J + span(extra)}."""
    Fp = core.reduce_mod(F, p); Tp = core.reduce_mod(T, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    nB = len(Bpolys)
    def space(mult, deg_tot, par_tot, extra_polys):
        rows = core.monomials(deg_tot, par_tot); idx = {m: k for k, m in enumerate(rows)}
        cols = []
        for i in range(6):
            pi = 0 if i < 3 else 1
            for m in core.monomials(deg_tot - 5, (par_tot - pi) % 2):
                cols.append(core.shift(Fi[i], m, 1, p))
        for e in extra_polys: cols.append(core.reduce_mod(e, p))
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
    E1, rkE = space(Tp, d + degT, (parD + parT_) % 2, [])
    Jr, rkJ = space(None, d, parD, extra or [])
    return E1, rkE, Jr, rkJ

def pivots(Mr, rk, n):
    out = []
    for i in range(rk):
        j = 0
        while j < n and int(Mr[i, j]) == 0: j += 1
        out.append(j)
    return out

def canon_complement(E1, rkE, Jr, rkJ, nB, p):
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

def inter_dim(Mr, rk, cols_keep, nB, p):
    if rk == 0: return 0
    comp = [j for j in range(nB) if j not in set(cols_keep)]
    if not comp: return rk
    sub_m = nmod_mat(len(comp), rk, [int(Mr[i, j]) for j in comp for i in range(rk)], p)
    _, nul = sub_m.nullspace()
    return nul

if __name__ == '__main__':
    t0 = time.time()
    bl = blocks(d)
    labels = []; Bpolys = []
    for lab, polys in bl:
        for pp_ in polys:
            labels.append(lab); Bpolys.append(pp_)
    nB = len(Bpolys)
    extra = []
    if d == 8:
        # lower kernel times linear forms (rational generators from degree 7 run)
        try:
            g7 = [core.parse(s) for s in json.load(open(f'{WORK}/kstruct2_{TARGET.replace("*","")}_d7.json'))['generators']]
            for G in g7:
                for i in range(3):
                    extra.append(core.mul(G, {tuple(1 if k == i else 0 for k in range(6)): Fraction(1)}))
            print('using', len(extra), 'products K_7 * X_i as extra trivial directions')
        except Exception as ex:
            print('no degree-7 generators available:', ex)
    print(f'TARGET={TARGET} d={d}: base-type space dim {nB}; blocks: {[(lab, len(pl)) for lab, pl in bl]}', flush=True)
    res = {}
    for p in (p1, p2):
        E1, rkE, Jr, rkJ = kernel_coords(Tq, Bpolys, p, extra)
        ker = canon_complement(E1, rkE, Jr, rkJ, nB, p)
        print(f'  p={p}: dim(E1 cap B)={rkE}, dim((J+extra) cap B)={rkJ}, kernel mod (J+extra) = {len(ker)}  ({time.time()-t0:.1f}s)', flush=True)
        res[p] = (E1, rkE, Jr, rkJ, ker)
        if p == p1:
            blab = [lab for lab, _ in bl]
            for r in range(1, len(blab) + 1):
                for sub in itertools.combinations(range(len(blab)), r):
                    keep = [j for j, lab in enumerate(labels) if lab in {blab[s] for s in sub}]
                    dE = inter_dim(E1, rkE, keep, nB, p); dJ = inter_dim(Jr, rkJ, keep, nB, p)
                    if dE - dJ > 0:
                        print(f'    blocks {[blab[s] for s in sub]}: E1-cap={dE}, J-cap={dJ}, kernel mod J = {dE - dJ}', flush=True)
    k1 = res[p1][4]; k2 = res[p2][4]
    gens = []; gens_struct = []
    for v1, v2 in zip(k1, k2):
        poly = {}; coords = []; ok = True
        for j in range(nB):
            a, mm = core.crt_pair(v1[j], p1, v2[j], p2)
            r = core.ratrecon(a, mm)
            if r is None: ok = False; break
            if r: coords.append((labels[j], j, r)); poly = core.add(poly, core.scale(Bpolys[j], r))
        if ok:
            print('  generator (structured coords):', [(lab, core.to_str({tuple(k for k in next(iter(Bpolys[j]))) : 1}) if False else lab, str(r)) for lab, j, r in coords][:60], '... nterms', len(coords), flush=True)
            gens.append(core.to_str(poly)); gens_struct.append([(lab, j, str(r)) for lab, j, r in coords])
        else:
            print('  ratrecon failed for a generator', flush=True)
    json.dump({'target': TARGET, 'd': d, 'labels': labels, 'generators': gens, 'structured': gens_struct,
               'basis': [core.to_str(b) for b in Bpolys]}, open(f'{WORK}/kstruct2_{TARGET.replace("*","")}_d{d}.json', 'w'))
    print('done', time.time() - t0)
