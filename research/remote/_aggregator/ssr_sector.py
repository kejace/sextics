"""Independent check: dimension of the S3xS3-invariant Y-odd sextics mod J(F) for symmetric-sextic-residual,
and of the (sgn,sgn)-polynomial sector (Delta_X Delta_Y). Residue character = det(g) * chi_A(g)."""
import sys, json, itertools, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core, jacobian
p = 2147483563
F, A, q, rec = core.load_instance('symmetric-sextic-residual', p)
t0=time.time()
mons = core.monomials(6, 1)  # Y-odd sextic monomials (exponent tuples)
def perm_apply(e, px, py):
    return (e[px[0]], e[px[1]], e[px[2]], 3+0 and e[3+py[0]], e[3+py[1]], e[3+py[2]])
def orbit_sum(e, sign=False):
    out = {}
    for px in itertools.permutations(range(3)):
        for py in itertools.permutations(range(3)):
            f = (e[px[0]], e[px[1]], e[px[2]], e[3+py[0]], e[3+py[1]], e[3+py[2]])
            c = 1
            if sign:
                sx = 1 if sum(1 for i in range(3) for j in range(i+1,3) if px[i]>px[j])%2==0 else -1
                sy = 1 if sum(1 for i in range(3) for j in range(i+1,3) if py[i]>py[j])%2==0 else -1
                c = sx*sy
            out[f] = (out.get(f,0) + c) % p
    return {k:v for k,v in out.items() if v}
seen=set(); inv=[]; sgn=[]
for e in mons:
    key = (tuple(sorted(e[:3])), tuple(sorted(e[3:])))
    if key in seen: continue
    seen.add(key)
    s = orbit_sum(e); 
    if s: inv.append(s)
    s2 = orbit_sum(e, sign=True)
    if s2: sgn.append(s2)
# reduce mod J and compute ranks of the spans of normal forms
from flint import nmod_mat
def rank_mod_J(polys):
    nfs = [jacobian.normal_form(F, P, p) for P in polys]
    nfs = [v for v in nfs if v]
    if not nfs: return 0, 0
    keys = sorted({k for v in nfs for k in v})
    idx = {k:i for i,k in enumerate(keys)}
    M = nmod_mat(len(nfs), len(keys), [0]*(len(nfs)*len(keys)), p)
    for r,v in enumerate(nfs):
        for k,c in v.items(): M[r, idx[k]] = c
    return M.rank(), len(nfs)
ri, ni = rank_mod_J(inv); rs, ns = rank_mod_J(sgn)
# A itself
nfA = jacobian.normal_form(F, A, p)
res = {'lane':'symmetric-sextic-residual','p':p,'invariant_Yodd_sextic_orbits':len(inv),'rank_mod_J_invariant':ri,
       'sgnsgn_Yodd_sextic_orbits':len(sgn),'rank_mod_J_sgnsgn':rs,'A_nonzero_mod_J': bool(nfA),
       'note':'residue character of an invariant polynomial = det = sgn(x)sgn(y); the target sector (residue sgn⊗sgn, tau=+1) has h^{3,1} = rank_mod_J_invariant; the polynomial-(sgn,sgn) sector (residue trivial) has h^{3,1} = rank_mod_J_sgnsgn',
       'seconds': round(time.time()-t0,1)}
print(json.dumps(res))
open('/tmp/hodge/results/_aggregator/spotcheck.jsonl','a').write(json.dumps(dict(res, kind='sector-check'))+'\n')
