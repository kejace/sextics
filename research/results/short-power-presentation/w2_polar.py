"""Polar residual spaces: PR_l(d) = (J_{d+1} ∩ l*S_d)/l ⊂ S_d for linear forms l, and their intersection with the E1
kernel K_d (mod J_d).  Usage: python3 w2_polar.py <d> <par> [--p P]"""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_families import form_polys, load_kernel_modp

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('par', type=int)
    ap.add_argument('--p', type=int, default=PGB[0]); a = ap.parse_args()
    F, A, q, rec = load(); d, par, p = a.d, a.par, a.p
    FP = form_polys(); FP['L+M'] = add(FP['L'], FP['M']); FP['L-M'] = sub(FP['L'], FP['M'])
    jref_d = JRef(F, d, par, p, verbose=False)
    sm = jref_d.std; sidx = {jref_d.mons[j]: t for t, j in enumerate(sm)}; ns = len(sm)
    Kp = load_kernel_modp(d, par, p) or []
    Knf = jref_d.nf_polys(Kp) if Kp else []
    def vec(nf):
        v = [0] * ns
        for m, c in nf.items(): v[sidx[m]] = c
        return v
    Kmat = nmod_mat(len(Knf), ns, [x for nf in Knf for x in vec(nf)], p) if Knf else None
    out = {}
    t0 = time.time()
    for name, lf in FP.items():
        lp = reduce_mod(lf, p); parl = core.parity_of(lp)
        pare = (par + parl) % 2   # parity of l*D
        gens = jacobian.jacobian_generators(F, d + 1, p, pare)
        Dm = monomials(d, par)
        lcols = [shift(lp, m, 1, p) for m in Dm]
        cols = [g[2] for g in gens] + lcols
        rows = {}
        for c in cols:
            for e in c: rows.setdefault(e, len(rows))
        R, C = len(rows), len(cols)
        flat = [0] * (R * C)
        for j, c in enumerate(cols):
            for e, v in c.items(): flat[rows[e] * C + j] = v
        M = nmod_mat(R, C, flat, p); del flat
        X, nul = M.nullspace()
        nJ = len(gens)
        Ds = []
        for c in range(nul):
            D = {}
            for j, m in enumerate(Dm):
                v = int(X[nJ + j, c])
                if v: D[m] = v
            if D: Ds.append(D)
        # mod J_d and intersection with K_d
        nfs = jref_d.nf_polys(Ds) if Ds else []
        nfs = [nf for nf in nfs if nf]
        rPR = nmod_mat(len(nfs), ns, [x for nf in nfs for x in vec(nf)], p).rank() if nfs else 0
        inter = None
        if nfs and Knf:
            both = nmod_mat(len(nfs) + len(Knf), ns, [x for nf in nfs for x in vec(nf)] + [x for nf in Knf for x in vec(nf)], p).rank()
            inter = rPR + len(Knf) - both
        out[name] = {'dim_PR_mod_J': rPR, 'dim_PR_incl_J': len(Ds), 'dim_PR_cap_K': inter}
        print(f'[polar] d={d} par={par} l={name}: dim (J_{d+1} ∩ l S_d)/l mod J_d = {rPR} (incl J: {len(Ds)}), ∩ K_d = {inter}', file=sys.stderr, flush=True)
    logrec = {'kind': 'polar-residuals', 'd': d, 'parity': par, 'p': p, 'dim_K_d': len(Knf), 'results': out, 'seconds': round(time.time() - t0, 1)}
    w2log(logrec); print(json.dumps(logrec), flush=True)

if __name__ == '__main__':
    main()
