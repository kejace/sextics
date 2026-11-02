"""dim Ann_R(l)_d for special linear forms l (mod p), d in a range, both parities; and whether Ann(l)_d ⊂ K_d (it must for l | A)."""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
from w2_families import form_polys
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dmin', type=int); ap.add_argument('dmax', type=int)
    ap.add_argument('--p', type=int, default=PGB[0]); ap.add_argument('--forms', default='X0,Y0,L,M,X1,Y1')
    a = ap.parse_args(); p = a.p
    F, A, q, rec = load(); FP = form_polys()
    out = {}
    t0 = time.time()
    for d in range(a.dmin, a.dmax + 1):
        for par in (0, 1):
            jref_d = JRef(F, d, par, p, verbose=False)
            for name in a.forms.split(','):
                lp = reduce_mod(FP[name], p); parl = core.parity_of(lp); pare = (par + parl) % 2
                gens = jacobian.jacobian_generators(F, d + 1, p, pare)
                Dm = monomials(d, par)
                cols = [g[2] for g in gens] + [shift(lp, m, 1, p) for m in Dm]
                rows = {}
                for c in cols:
                    for e in c: rows.setdefault(e, len(rows))
                R, C = len(rows), len(cols)
                flat = [0] * (R * C)
                for j, c in enumerate(cols):
                    for e, v in c.items(): flat[rows[e] * C + j] = v
                M = nmod_mat(R, C, flat, p); del flat
                X, nul = M.nullspace(); nJ = len(gens)
                Ds = []
                for c in range(nul):
                    D = {Dm[j]: int(X[nJ + j, c]) for j in range(len(Dm)) if int(X[nJ + j, c])}
                    if D: Ds.append(D)
                nfs = [nf for nf in (jref_d.nf_polys(Ds) if Ds else []) if nf]
                sm = jref_d.std; sidx = {jref_d.mons[j]: t for t, j in enumerate(sm)}
                rk = nmod_mat(len(nfs), len(sm), [nf.get(jref_d.mons[sm[t]], 0) for nf in nfs for t in range(len(sm))], p).rank() if nfs else 0
                out[f'd{d}p{par}:{name}'] = rk
                print(f'[annl] d={d} par={par} l={name}: dim Ann_R(l)_d = {rk}  ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    logrec = {'kind': 'ann-linear-forms', 'p': p, 'dims': out, 'seconds': round(time.time() - t0, 1)}
    w2log(logrec); print(json.dumps(logrec), flush=True)
if __name__ == '__main__':
    main()
