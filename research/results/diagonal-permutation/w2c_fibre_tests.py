"""J1'/J2': hidden structure on the (01)-pencil fibre X_mu = X ∩ {X0-X1 = mu (Y0-Y1)} for a random rational mu.
Coordinates on the hyperplane: (p,q,a,b,t) with X0=(p+mu t)/2, X1=(p-mu t)/2, X2=a, Y0=(q+t)/2, Y1=(q-t)/2, Y2=b.
F_mu is even in t (the reflection s).  J1': points v in P^4 with rank Hess F_mu(v) <= 2 (needed for a §2.7 exchange partner
coordinate (rank 1, v off X_mu) or a cone tangent section).  J2': hidden reflections of F_mu (cubic polar condition), i.e. a
second reflection on the fibre, which would make the fibrewise quotient Fano.  Both via msolve over F_p, all affine charts."""
import sys, json, time, subprocess, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import sympy as sp
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
p_,q_,a_,b_,t_ = sp.symbols('p q a b t'); W = [p_,q_,a_,b_,t_]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.sympify(core.to_str(F))
mu = sp.Rational(3, 7)
Fmu = sp.expand(Fs.subs({X0:(p_+mu*t_)/2, X1:(p_-mu*t_)/2, X2:a_, Y0:(q_+t_)/2, Y1:(q_-t_)/2, Y2:b_}, simultaneous=True))
assert sp.expand(Fmu.subs(t_, -t_) - Fmu) == 0
pr = 1073741827
nthreads = int(os.environ.get('NTHREADS','16'))
def to_ms(poly, unknowns):
    poly = sp.Poly(poly, *unknowns)
    den = 1
    for c in poly.coeffs(): den = sp.ilcm(den, sp.Rational(c).q)
    poly = sp.Poly(sp.expand(poly.as_expr()*den), *unknowns)
    terms = []
    for mon, c in poly.terms():
        c = int(c) % pr
        if c == 0: continue
        mstr = '*'.join('%s^%d'%(str(unknowns[i]),e) for i,e in enumerate(mon) if e>0)
        terms.append(('%d*%s'%(c, mstr)) if mstr else '%d'%c)
    return '+'.join(terms) if terms else '0'
def run(fname, unknowns, eqs, tag):
    with open(fname,'w') as f:
        f.write(','.join(str(x) for x in unknowns)+'\n'); f.write('%d\n'%pr)
        f.write(',\n'.join(to_ms(e, unknowns) for e in eqs if sp.expand(e)!=0)+'\n')
    t1=time.time()
    r = subprocess.run(['msolve','-v','1','-t',str(nthreads),'-f',fname,'-o',fname.replace('.ms','.out')], capture_output=True, text=True, timeout=7000)
    out = open(fname.replace('.ms','.out')).read().strip()
    res = {"exp":tag,"file":fname,"prime":pr,"mu":str(mu),"result":("empty (GB={1})" if out.startswith('[1') or out.strip()=='[]' else "see "+fname.replace('.ms','.out')),
           "raw_head":out[:300],"seconds":round(time.time()-t1,1),"job":"construct-dp-fibre","stdout_tail":r.stdout[-400:]}
    print(res, flush=True)
    with open('wave2-construct.log.jsonl','a') as f: f.write(json.dumps(res)+'\n')
which = sys.argv[1]
if which == 'J1f':
    H = sp.Matrix(5,5, lambda i,j: sp.expand(sp.diff(Fmu, W[i], W[j])))
    minors = []
    for rows in itertools.combinations(range(5),3):
        for cols in itertools.combinations(range(5),3):
            if rows <= cols:
                m = sp.expand(H.extract(list(rows), list(cols)).det())
                if m != 0: minors.append(m)
    for ch in range(5):
        eqs = [sp.expand(m.subs(W[ch], 1)) for m in minors]
        unknowns = [w for i,w in enumerate(W) if i != ch]
        run('fibre_hessrank_chart%d.ms'%ch, unknowns, eqs, 'J1f-fibre-hessrank2 chart %s=1'%W[ch])
elif which == 'J2f':
    v = sp.symbols('v0:5'); u = sp.Symbol('u'); qs = sp.symbols('q0:15')
    quad = [W[i]*W[j] for i in range(5) for j in range(i,5)]
    for ch in range(5):
        vv = list(v); vv[ch] = 1
        def dv(e): return sp.expand(sum(vv[i]*sp.diff(e, W[i]) for i in range(5)))
        d3 = dv(dv(dv(Fmu)))
        lv = sp.expand(sum(sp.diff(Fmu, W[i]).subs(dict(zip(W, vv)))*W[i] for i in range(5)))
        Fv = sp.expand(Fmu.subs(dict(zip(W, vv))))
        Q = sum(qs[k]*quad[k] for k in range(15))
        E = sp.Poly(sp.expand(d3 - lv*Q), *W)
        eqs = list(E.coeffs()) + [sp.expand(u*Fv - 1)]
        unknowns = [v[i] for i in range(5) if i != ch] + list(qs) + [u]
        run('fibre_reflection_chart%d.ms'%ch, unknowns, eqs, 'J2f-fibre-hidden-reflection chart v%d=1'%ch)
