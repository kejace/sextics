"""C1: bigrading criterion for a sect.2.7 'polar cover with exchange'.
A cover z <-> T, T^6 = P(w'') (w'' a block of coordinates, X = rest) makes F symmetric in (z,T) iff,
in coordinates (z, w'', X), every monomial of F has deg_z = deg_{w''} except z^6 and the pure w''-sextic P.
(Proof: F = sum z^k T^m g_{km}(X, base); g_{km}=g_{mk} and polynomiality force F_{k,m}=0 unless k=m or {k,m}={0,6}.)
We test every coordinate system in a list, every coordinate z, every block w'' of the others (size 1..4).
"""
import sys, json, itertools, time
import sympy as sp
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
F, A, q, rec = core.load_instance('diagonal-permutation')
X0,X1,X2,Y0,Y1,Y2 = Z = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
Fs = sp.sympify(core.to_str(F).replace('^','**'), locals=dict(zip(['X0','X1','X2','Y0','Y1','Y2'], Z)))
w = sp.exp(2*sp.pi*sp.I/3)
# coordinate systems: dict name -> (new symbols, substitution Z_i = expr(new))
p,r,a,q_,s,b = sp.symbols('p r a q s b')
P,R1,R2,Q,S1,S2 = sp.symbols('P R1 R2 Q S1 S2')
xi0,xi1,xi2,eta0,eta1,eta2 = sp.symbols('xi0 xi1 xi2 eta0 eta1 eta2')
systems = {
 'coord': (Z, dict(zip(Z, Z))),
 '01-adapted (p,r,a,q,s,b)': ((p,r,a,q_,s,b), {X0:(p+r)/2, X1:(p-r)/2, X2:a, Y0:(q_+s)/2, Y1:(q_-s)/2, Y2:b}),
 'S3-adapted (P,R1,R2,Q,S1,S2)': ((P,R1,R2,Q,S1,S2), {X0:P+R1, X1:P+R2, X2:P-R1-R2, Y0:Q+S1, Y1:Q+S2, Y2:Q-S1-S2}),
 'A3-Fourier (xi,eta)': ((xi0,xi1,xi2,eta0,eta1,eta2), {X0:(xi0+xi1+xi2)/3, X1:(xi0+w**2*xi1+w*xi2)/3, X2:(xi0+w*xi1+w**2*xi2)/3,
                                               Y0:(eta0+eta1+eta2)/3, Y1:(eta0+w**2*eta1+w*eta2)/3, Y2:(eta0+w*eta1+w**2*eta2)/3}),
}
out = []
t0 = time.time()
for name,(V,sub) in systems.items():
    G = sp.Poly(sp.expand(Fs.subs(sub, simultaneous=True)), *V)
    G = sp.Poly(sp.nsimplify(sp.expand(G.as_expr()), rational=True) if name!='A3-Fourier (xi,eta)' else sp.expand(G.as_expr()), *V)
    terms = G.terms()
    nterms = len(terms)
    best = []
    for zi in range(6):
        others = [j for j in range(6) if j!=zi]
        for k in range(1,5):
            for blk in itertools.combinations(others, k):
                viol = 0; viol_list=[]
                for mono, c in terms:
                    dz = mono[zi]; dw = sum(mono[j] for j in blk)
                    if dz == dw: continue
                    if dz == 6 and dw == 0: continue
                    if dz == 0 and dw == 6: continue
                    viol += 1
                    if len(viol_list) < 4: viol_list.append(''.join(f'{V[j]}^{mono[j]}' for j in range(6) if mono[j]))
                best.append((viol, str(V[zi]), [str(V[j]) for j in blk], viol_list))
    best.sort(key=lambda t: t[0])
    out.append({'system': name, 'nterms': nterms, 'nconfigs': len(best), 'min_violations': best[0][0], 'best5': best[:5]})
    print(name, 'terms', nterms, 'configs', len(best), 'min viol', best[0][0], 'best', best[:3], flush=True)
res = {'exp': 'C1-polar-cover-bigrading', 'criterion': 'deg_z == deg_w'' for all monomials except z^6 and pure w''-sextic',
       'results': out, 'seconds': round(time.time()-t0,1), 'verdict': 'no (z, block) in any tested coordinate system satisfies the criterion' if all(o['min_violations']>0 for o in out) else 'FOUND'}
json.dump(res, open('c1_polar_cover.json','w'), indent=1)
with open('../wave2-construct.log.jsonl','a') as f: f.write(json.dumps({k:v for k,v in res.items() if k!='results'} | {'min_violations': {o['system']:o['min_violations'] for o in out}, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')})+'\n')
print('DONE', res['verdict'])
