#!/usr/bin/env python3
"""Generate msolve inputs for three projective loci attached to the dense-rational sextic G = 9973 F:
  V1 = {e : rank Hess G(e) <= 1}      (2x2 minors, degree 8)   -- holdout-type exchange covers, order-3/6 quasi-reflections
  V2 = {e : rank Hess G(e) <= 2}      (3x3 minors, degree 12)  -- all quasi-reflections (any order), all exchange covers
  V3 = {e : D^3G(e)|_{ker dG(e)} = 0} (degree 15)              -- linear reflections (G even in some linear form)
      T'_{ijk} = 27 G^2 G_{ijk} - 18 G (G_i G_{jk} + G_j G_{ik} + G_k G_{ij}) + 30 G_i G_j G_k  (i<=j<=k).
Emptiness mod p of a homogeneous system (pure power of every variable in the DRL leading ideal) is rigorous over C.
"""
import sys, json, time, itertools, subprocess
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint
P = 1073741789
OUT = '/tmp/hodge/results/dense-rational/wave2-construct'
LOG = '/tmp/hodge/results/dense-rational/wave2-construct.log.jsonl'
def log(**kw):
    kw['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    with open(LOG, 'a') as f: f.write(json.dumps(kw) + '\n')
F, A, q, rec = core.load_instance('dense-rational')
ctx = flint.fmpz_mpoly_ctx.get(('e0','e1','e2','e3','e4','e5'), 'degrevlex')
E = ctx.gens()
def mono(ex):
    m = ctx.from_dict({tuple(ex): 1}); return m
G = ctx.from_dict({tuple(ex): int(c * 9973) for ex, c in F.items()})
assert all((c * 9973).denominator == 1 for c in F.values())
import os
PFX = os.environ.get('PFX', '')
if os.environ.get('FERMAT'):
    G = sum(E[i]**6 for i in range(6))
t0 = time.time()
Gi = [G.derivative(i) for i in range(6)]
H = [[Gi[i].derivative(j) for j in range(6)] for i in range(6)]
def fmt(poly):
    d = poly.to_dict()
    terms = []
    for ex, c in d.items():
        c = int(c) % P
        if c == 0: continue
        mon = '*'.join(f'e{i}^{k}' if k > 1 else f'e{i}' for i, k in enumerate(ex) if k > 0)
        terms.append(f'{c}*{mon}' if mon else f'{c}')
    return '+'.join(terms)
def write(name, polys, note):
    seen, eqs = set(), []
    for pz in polys:
        if pz.is_zero(): continue
        s = fmt(pz)
        if s and s not in seen: seen.add(s); eqs.append(s)
    fn = f'{OUT}/{PFX}{name}_p{P}.ms'
    with open(fn, 'w') as f:
        f.write('e0,e1,e2,e3,e4,e5\n'); f.write(f'{P}\n'); f.write(',\n'.join(eqs) + '\n')
    degs = sorted(set(int(pz.total_degree()) for pz in polys if not pz.is_zero()))
    log(exp='msolve-input', locus=name, file=fn, prime=P, n_equations=len(eqs), degrees=degs, note=note, seconds=round(time.time()-t0,1))
    print(name, len(eqs), 'equations, degrees', degs, flush=True)
    return fn
which = sys.argv[1:] or ['V1', 'V2', 'V3']
if 'V1' in which:
    polys = [H[i][j]*H[k][l] - H[i][l]*H[k][j] for i, k in itertools.combinations(range(6), 2) for j, l in itertools.combinations(range(6), 2)]
    write('V1_hessrank1', polys, '2x2 minors of Hess G: rank<=1 locus (exchange covers / order-3,6 quasi-reflections)')
if 'V2' in which:
    def det3(r, c):
        (a,b,cc),(d,e,f) = r, c
        M = [[H[x][y] for y in c] for x in r]
        return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    polys = [det3(r, c) for r in itertools.combinations(range(6), 3) for c in itertools.combinations(range(6), 3)]
    write('V2_hessrank2', polys, '3x3 minors of Hess G: rank<=2 locus (all quasi-reflections, all holdout-type exchange covers)')
if 'V3' in which:
    G2 = G*G
    polys = []
    for i, j, k in itertools.combinations_with_replacement(range(6), 3):
        Gijk = H[i][j].derivative(k)
        T = 27*G2*Gijk - 18*G*(Gi[i]*H[j][k] + Gi[j]*H[i][k] + Gi[k]*H[i][j]) + 30*Gi[i]*Gi[j]*Gi[k]
        polys.append(T)
    write('V3_reflection', polys, 'cubic Taylor term at e vanishes on the polar hyperplane: necessary for G even in a linear form (linear reflection)')
print('done', round(time.time()-t0,1), 's')
