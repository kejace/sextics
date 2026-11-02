"""dim Ann_R(m)_d for monomial multipliers m dividing A = X0^3 Y0^3 (mod p). Usage: python3 w2_annpow.py dmin dmax [--maxk K]"""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dmin', type=int); ap.add_argument('dmax', type=int)
    ap.add_argument('--maxk', type=int, default=5); ap.add_argument('--p', type=int, default=PGB[0]); a = ap.parse_args(); p = a.p
    F, A, q, rec = load()
    mults = []
    for i in range(4):
        for j in range(4):
            if 0 < i + j <= a.maxk: mults.append((i, j))
    out = {}; t0 = time.time()
    for d in range(a.dmin, a.dmax + 1):
        for par in (0, 1):
            jref_d = JRef(F, d, par, p, verbose=False)
            sm = jref_d.std; sidx = {jref_d.mons[j]: t for t, j in enumerate(sm)}
            jrefs = {}
            for (i, j) in mults:
                mexp = (i, 0, 0, j, 0, 0); k = i + j; parm = j % 2; pare = (par + parm) % 2
                if (d + k, pare) not in jrefs:
                    jrefs[(d + k, pare)] = JRef(F, d + k, pare, p, verbose=False)
                jr = jrefs[(d + k, pare)]
                Dm = monomials(d, par)
                # kernel of multiplication by the monomial: NF(m * mon) for all mon -> left nullspace (same trick as e1_kernel)
                mm = [mon_add(m_, mexp) for m_ in Dm]
                Pk = [t for t, mu in enumerate(mm) if jr.idx[mu] in jr.pivset]; Tk = [t for t, mu in enumerate(mm) if jr.idx[mu] not in jr.pivset]
                StdA = {jr.idx[mm[t]]: t for t in Tk}; std_other = [c for c in jr.std if c not in StdA]
                rows = [jr.row(jr.rowof[jr.idx[mm[t]]]) for t in Pk]
                nP = len(Pk)
                if nP == 0:
                    dimk = 0
                else:
                    flat = [0] * (len(std_other) * nP)
                    for a_, r in enumerate(rows):
                        for t, c in enumerate(std_other):
                            if r[c]: flat[t * nP + a_] = r[c]
                    E = nmod_mat(len(std_other), nP, flat, p); X, nul = E.nullspace()
                    vecs = []
                    for c in range(nul):
                        xP = [int(X[a_, c]) for a_ in range(nP)]
                        D = {Dm[Pk[a_]]: xP[a_] for a_ in range(nP) if xP[a_]}
                        for cidx, t in StdA.items():
                            v = sum(xP[a_] * rows[a_][cidx] for a_ in range(nP)) % p
                            if v: D[Dm[t]] = v
                        vecs.append(D)
                    nfs = [nf for nf in (jref_d.nf_polys(vecs) if vecs else []) if nf]
                    dimk = nmod_mat(len(nfs), len(sm), [nf.get(jref_d.mons[sm[t]], 0) for nf in nfs for t in range(len(sm))], p).rank() if nfs else 0
                out[f'd{d}p{par}:X0^{i}Y0^{j}'] = dimk
                print(f'[annpow] d={d} par={par} m=X0^{i}Y0^{j}: dim Ann_R(m)_d mod J = {dimk} ({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
    logrec = {'kind': 'ann-monomials', 'p': p, 'dims': out, 'seconds': round(time.time() - t0, 1)}
    w2log(logrec); print(json.dumps(logrec), flush=True)
if __name__ == '__main__':
    main()
