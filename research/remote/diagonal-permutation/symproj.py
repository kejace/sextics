"""S3 / transposition projections of a linear family of polynomials (exact)."""
import sys; sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
from hodge.core import add, scale
from fractions import Fraction
from flint import fmpq_mat, fmpq
PERMS = [(0,1,2),(1,0,2),(0,2,1),(2,1,0),(1,2,0),(2,0,1)]
def act(P, s):
    """(s.P)(Z) = P(Z_{s}) : permute indices of X and Y simultaneously."""
    out = {}
    for e, c in P.items():
        e2 = (e[s[0]], e[s[1]], e[s[2]], e[3+s[0]], e[3+s[1]], e[3+s[2]])
        out[e2] = out.get(e2, 0) + c
    return {e: v for e, v in out.items() if v}
def average(P, perms=PERMS, signs=None):
    out = {}
    for k, s in enumerate(perms):
        out = add(out, scale(act(P, s), Fraction(1 if signs is None else signs[k], len(perms))))
    return out
def span_basis(polys):
    polys = [P for P in polys if P]
    if not polys: return []
    allm = {}
    for P in polys:
        for e in P: allm.setdefault(e, len(allm))
    M = fmpq_mat(len(polys), len(allm))
    for r, P in enumerate(polys):
        for e, v in P.items(): M[r, allm[e]] = fmpq(v.numerator, v.denominator)
    Mr, rk = M.rref()
    inv = {v: e for e, v in allm.items()}
    return [{inv[c]: Fraction(int(Mr[i, c].p), int(Mr[i, c].q)) for c in range(len(allm)) if Mr[i, c] != 0} for i in range(rk)]
def sign_of(s):
    # parity of permutation
    inv = sum(1 for i in range(3) for j in range(i+1,3) if s[i] > s[j])
    return -1 if inv % 2 else 1
def invariant_part(polys): return span_basis([average(P) for P in polys])
def sign_part(polys): return span_basis([average(P, signs=[sign_of(s) for s in PERMS]) for P in polys])
def t01_invariant(polys): return span_basis([average(P, perms=[(0,1,2),(1,0,2)]) for P in polys])
def t01_anti(polys): return span_basis([average(P, perms=[(0,1,2),(1,0,2)], signs=[1,-1]) for P in polys])
if __name__ == '__main__':
    import json
    fam = json.load(open('/tmp/hodge/results/diagonal-permutation/families.json'))
    out = {}
    for key, strs in fam.items():
        polys = [core.parse(s) for s in strs]
        inv = invariant_part(polys); sg = sign_part(polys); ti = t01_invariant(polys); ta = t01_anti(polys)
        print(f'{key}: dim={len(polys)} S3-inv={len(inv)} S3-sign={len(sg)} (01)-inv={len(ti)} (01)-anti={len(ta)}')
        out[key] = {'inv': [core.to_str(P) for P in inv], 'sign': [core.to_str(P) for P in sg],
                    't01inv': [core.to_str(P) for P in ti], 't01anti': [core.to_str(P) for P in ta]}
    json.dump(out, open('/tmp/hodge/results/diagonal-permutation/families_sym.json', 'w'), indent=1)
