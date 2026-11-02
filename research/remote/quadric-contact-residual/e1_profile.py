"""E1 kernel dims per G_aff character at degree d: dim{D in S_d[c]: A D in J} - dim J_d[c]."""
import sys, json, time
from inv import *
d = int(sys.argv[1])
t0 = time.time()
chars = sorted({chi(e) for e in monomials(d)})
out = {}
tot_extra = 0
extra_dirs = {}
for c in chars:
    mons = [e for e in monomials(d) if chi(e) == c]
    S6 = Sector(d + 6, chi_add(c, CHI_A))
    vecs = [S6.vec(shift(Ap, e, 1, p)) for e in mons]
    NF = S6.nf_matrix(vecs)
    X, nul = NF.transpose().nullspace()
    # dim J_d[c]
    Sd = Sector(d, c)
    extra = nul - Sd.rank
    out[str(c)] = {'dimS': len(mons), 'dimJ': Sd.rank, 'kernel': nul, 'extra': extra}
    tot_extra += extra
    if extra > 0:
        # extract kernel vectors modulo J_d[c]: reduce kernel polys with Sd normal form, keep independent
        kers = [[int(X[i, j]) for i in range(len(mons))] for j in range(nul)]
        polys = []
        for v in kers:
            P = {e: cv for e, cv in zip(mons, v) if cv}
            polys.append(P)
        NFk = Sd.nf_matrix([Sd.vec(P) for P in polys])
        R, rk = NFk.rref()
        # pick independent kernel polys
        chosen = []
        cur = []
        for i, P in enumerate(polys):
            test = nmod_mat(len(cur) + 1, NFk.ncols(), [int(NFk[t, j]) for t in cur + [i] for j in range(NFk.ncols())], p)
            if test.rank() == len(cur) + 1:
                cur.append(i); chosen.append(P)
            if len(cur) == extra:
                break
        extra_dirs[str(c)] = [core.to_str(P) for P in chosen]
    print(c, out[str(c)], f'({time.time()-t0:.0f}s)', file=sys.stderr, flush=True)
summary = {'d': d, 'total_extra': tot_extra, 'dimJ_total': sum(v['dimJ'] for v in out.values()), 'per_char': out, 'seconds': round(time.time() - t0, 1)}
json.dump({'summary': summary, 'extra_dirs': extra_dirs}, open(f'e1_profile_d{d}.json', 'w'))
print(json.dumps({k: v for k, v in summary.items() if k != 'per_char'}))
for c, L in extra_dirs.items():
    print('extra', c, len(L))
