"""Ideal of the eight planes in degree 10: (Q)*S_8 and (Q)^2*S_6 (both parities) vs the admissible space (mod p)."""
import sys, json, time, pickle, itertools
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *
import families as fam_mod
from hodge.core import monomials, mul, add, reduce_mod, parity_of
from flint import nmod_mat
p = core.PRIMES[1]
d = 10
def rank_of(polys, monos):
    if not polys: return 0
    R, C = len(monos), len(polys); flat = [0] * (R * C)
    for j, P_ in enumerate(polys):
        for e, v in P_.items(): flat[monos[e] * C + j] = v % p
    return nmod_mat(R, C, flat, p).rank()
Qp = [reduce_mod(Q0, p), reduce_mod(Q1, p), reduce_mod(Q2, p)]
out = []
for par in (0, 1):
    kernel = pickle.load(open(f'{WD}/e1kernel_basis_d{d}_par{par}.pkl', 'rb'))
    J = fam_mod.J_basis(d, par, p); adm = kernel + J
    fams = {}
    fams['(Q)*S_8'] = [mul(Qp[j], {m: 1}, p) for j in range(3) for m in monomials(8, par)]
    fams['(Q)^2*S_6'] = [mul(mul(Qp[j], Qp[k], p), {m: 1}, p) for j in range(3) for k in range(j, 3) for m in monomials(6, par)]
    monos = {}
    for lst in [adm] + list(fams.values()):
        for P_ in lst:
            for e in P_: monos.setdefault(e, len(monos))
    rA = rank_of(adm, monos); rJ = rank_of(J, monos)
    for name, lst in fams.items():
        t0 = time.time()
        rF = rank_of(lst, monos); rAF = rank_of(adm + lst, monos); rJF = rank_of(J + lst, monos)
        rec = {'d': d, 'par': par, 'family': name, 'n_gens': len(lst), 'dim_fam': rF, 'dim_S': len(monomials(d, par)), 'fam_cap_admissible': rF + rA - rAF,
               'fam_cap_J': rF + rJ - rJF, 'new_directions': (rF + rA - rAF) - (rF + rJ - rJF), 'dim_adm': rA, 'dim_J': rJ, 'seconds': round(time.time() - t0, 1)}
        print(json.dumps(rec), flush=True); out.append(rec)
        with open(f'{WD}/wave2-search.log.jsonl', 'a') as f:
            f.write(json.dumps({'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'eight-plane ideal vs admissible (local, mod p)', **rec}) + '\n')
json.dump(out, open(f'{WD}/fam10b-results.json', 'w'))
