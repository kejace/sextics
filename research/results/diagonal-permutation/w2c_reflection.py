"""Wave-2 construct lens, job J2: hidden reflections of F.
A linear involution r = reflection (eigenvalues (1^5,-1)) with F o r = F has -1-eigenvector v with F(v) != 0,
fixed hyperplane = polar hyperplane ker(l_v), l_v(z) = sum_i dF/dz_i(v) z_i, and the odd polars vanish on it:
l_v | d_v^3 F (cubic) and l_v | d_v F (quintic).  We solve the cubic condition with a cofactor quadric q:
d_v^3 F = l_v * q  (56 coefficient equations, unknowns v (chart) and q (21)), plus u*F(v)-1=0.
"""
import sys, json, time, subprocess, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import sympy as sp
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
Z = [X0,X1,X2,Y0,Y1,Y2]
v = sp.symbols('v0:6'); u = sp.Symbol('u')
qs = sp.symbols('q0:21')
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.sympify(core.to_str(F))
p = 1073741827
chart = int(sys.argv[1]) if len(sys.argv)>1 else 0
onX = (len(sys.argv)>2 and sys.argv[2]=='onX')   # variant: v on X (cone tangent sections) instead of F(v)!=0
vv = list(v); vv[chart] = 1
# directional derivative operator along v
def dv(e):
    return sp.expand(sum(vv[i]*sp.diff(e, Z[i]) for i in range(6)))
d1 = dv(Fs); d3 = dv(dv(d1))
lv = sp.expand(sum(sp.diff(Fs, Z[i]).subs(dict(zip(Z, vv)))*Z[i] for i in range(6)))
Fv = sp.expand(Fs.subs(dict(zip(Z, vv))))
quad_monos = [Z[i]*Z[j] for i in range(6) for j in range(i,6)]
Q = sum(qs[k]*quad_monos[k] for k in range(21))
E = sp.Poly(sp.expand(d3 - lv*Q), *Z)
eqs = [c for c in E.coeffs()]
if onX:
    eqs.append(sp.expand(Fv))
    unknowns = [v[i] for i in range(6) if i != chart] + list(qs)
else:
    eqs.append(sp.expand(u*Fv - 1))
    unknowns = [v[i] for i in range(6) if i != chart] + list(qs) + [u]
def to_ms(poly):
    poly = sp.Poly(poly, *unknowns)
    den = 1
    for c in poly.coeffs():
        den = sp.ilcm(den, sp.Rational(c).q)
    poly = sp.Poly(sp.expand(poly.as_expr()*den), *unknowns)
    terms = []
    for mon, c in poly.terms():
        c = int(c) % p
        if c == 0: continue
        mstr = '*'.join('%s^%d'%(str(unknowns[i]),e) for i,e in enumerate(mon) if e>0)
        terms.append(('%d*%s'%(c, mstr)) if mstr else '%d'%c)
    return '+'.join(terms) if terms else '0'
tag = ('onX' if onX else 'offX'); fname = 'reflection_%s_v%d.ms'%(tag,chart)
with open(fname,'w') as f:
    f.write(','.join(str(x) for x in unknowns)+'\n'); f.write('%d\n'%p)
    f.write(',\n'.join(to_ms(e) for e in eqs if sp.expand(e)!=0)+'\n')
print('wrote', fname, len(eqs), 'equations', len(unknowns), 'unknowns', flush=True)
t1=time.time()
nthreads = int(os.environ.get('NTHREADS','16'))
r = subprocess.run(['msolve','-v','1','-t',str(nthreads),'-f',fname,'-o',fname.replace('.ms','.out')], capture_output=True, text=True, timeout=7000)
print(r.stdout[-3000:], r.stderr[-3000:], flush=True)
out = open(fname.replace('.ms','.out')).read().strip()
print('OUTPUT:', out[:2000], flush=True)
res = {"exp":"J2-cubic-polar-condition","family":"points v with l_v | d_v^3 F, %s (chart v%d=1)"%(("v on X: cone tangent sections" if onX else "F(v)!=0: hidden reflections / §2.7 partner"),chart),"prime":p,
       "result": ("none (GB={1})" if out.startswith('[1') or out.strip()=='[]' else "see "+fname.replace(".ms",".out")),
       "raw_head": out[:200], "seconds": round(time.time()-t1,1), "job":"construct-dp-j1j2",
       "note":"cubic-polar condition l_v | d_v^3 F with F(v)!=0; empty => no hidden reflection (mod p evidence)"}
with open('wave2-construct.log.jsonl','a') as f: f.write(json.dumps(res)+'\n')
print(res)
