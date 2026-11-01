"""target_image analysis: which degree-6 Y-odd targets does D serve? Decompose by X2-exponent mod 3 sector and
by (X-degree, Y-degree) bidegree; test whether A is in the served space (mod p); report dim modulo J_6."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import ctr
from ctr import *
from flint import nmod_mat
p = core.PRIMES[0]

def served(D, F_=None, label=''):
    F_ = F if F_ is None else F_
    t0 = time.time()
    dim, basis, S = fixed_d.target_image(F_, 1, D, p=p)
    mons = core.monomials(6, 1)
    idx = {m: i for i, m in enumerate(mons)}
    # J_6^odd dimension and whether image contains J_6
    Fp = core.reduce_mod(F_, p)
    Fi = [core.pderiv(Fp, i, p) for i in range(6)]
    Jcols = []
    for i in range(6):
        par = 1 if i < 3 else 0
        for m in core.monomials(1, par):
            Jcols.append(core.shift(Fi[i], m, 1, p))
    def mat(polys):
        flat = [0] * (len(polys) * len(mons))
        for j, P_ in enumerate(polys):
            for e, v in P_.items(): flat[j * len(mons) + idx[e]] = v
        return nmod_mat(len(polys), len(mons), flat, p)
    MJ = mat(Jcols); rJ = MJ.rank()
    MB = mat(basis) if basis else None
    rB = MB.rank() if basis else 0
    both = mat(basis + Jcols) if basis else MJ
    rBoth = both.rank()
    Ap = core.reduce_mod(A, p)
    rA = mat(basis + [Ap]).rank() if basis else 1
    rAJ = mat(basis + Jcols + [Ap]).rank()
    out = {'label': label, 'dim_image': rB, 'dim_J6odd': rJ, 'dim_image+J': rBoth, 'image_mod_J': rBoth - rJ,
           'A_in_image': rA == rB, 'A_in_image+J': rAJ == rBoth, 'seconds': round(time.time() - t0, 1)}
    # sector decomposition of the image modulo J: project basis onto X2-exponent mod 3 sectors is not linear-invariant
    # unless D is g-invariant; report the sector content of a basis anyway (support statistics)
    sect = [0, 0, 0]
    bid = {}
    for b in basis:
        for e in b:
            sect[e[2] % 3] += 1
            k = (e[0]+e[1]+e[2], e[3]+e[4]+e[5]); bid[k] = bid.get(k, 0) + 1
    out['basis_support_by_X2mod3'] = sect
    out['basis_support_by_bidegree'] = {str(k): v for k, v in sorted(bid.items())}
    return out, basis

if __name__ == '__main__':
    for Dexpr in sys.argv[1:]:
        D = P(Dexpr)
        out, basis = served(D, label=Dexpr)
        print(json.dumps(out), flush=True)
        log({'engine': 2, 'D': Dexpr[:400], 'd': core.degree(D), 'family': 'target_image', 'e1': None, 'consistent': None,
             'rank': out['dim_image'], 'nullity': out['image_mod_J'], 'seconds': out['seconds'],
             'note': f"served space dim {out['dim_image']} (mod J: {out['image_mod_J']}); A served: {out['A_in_image+J']}"})
