"""Wave-2 construct lens, job J1: locus of points v in P^5 where Hess(F)(v) has rank <= 2.
rank 1 at v not in X  <=> F = c y^6 + y G5 + G6 in coordinates with y dual to v  (necessary for a §2.7-type
exchange cover with partner coordinate y);  rank <= 2 at v in X <=> tangent hyperplane section is a cone with vertex v
(necessary).  We compute a Groebner basis of the ideal of all 3x3 minors over F_p with msolve.
Output: hessrank_<chart>.out ; summary line in wave2-construct.log.jsonl
"""
import sys, json, time, subprocess, itertools, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import sympy as sp
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
V = [X0,X1,X2,Y0,Y1,Y2]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.sympify(core.to_str(F))
p = 1073741827
H = sp.Matrix(6,6, lambda i,j: sp.expand(sp.diff(Fs, V[i], V[j])))
t0=time.time()
minors = []
for rows in itertools.combinations(range(6),3):
    for cols in itertools.combinations(range(6),3):
        if rows <= cols:  # symmetric matrix: minor(rows,cols) = minor(cols,rows)
            m = sp.expand(H.extract(list(rows), list(cols)).det())
            if m != 0:
                minors.append(m)
print('minors', len(minors), 'time', time.time()-t0, flush=True)
chart = sys.argv[1] if len(sys.argv) > 1 else 'X0'
charts = {'X0':X0,'X1':X1,'X2':X2,'Y0':Y0,'Y1':Y1,'Y2':Y2}
def to_ms(poly):
    # msolve wants integer coefficients; clear denominators
    poly = sp.Poly(poly, *V)
    den = 1
    for c in poly.coeffs():
        den = sp.ilcm(den, sp.Rational(c).q)
    poly = sp.Poly(sp.expand(poly.as_expr()*den), *V)
    terms = []
    for mon, c in poly.terms():
        c = int(c) % p
        if c == 0: continue
        mstr = '*'.join('%s^%d'%(str(V[i]),e) for i,e in enumerate(mon) if e>0)
        terms.append(('%d*%s'%(c, mstr)) if mstr else '%d'%c)
    return '+'.join(terms) if terms else '0'
eqs = [to_ms(m) for m in minors]
# chart: set chosen variable = 1 (dehomogenize) -> substitute
eqs_chart = []
for m in minors:
    mm = sp.expand(m.subs(charts[chart], 1))
    if mm != 0:
        eqs_chart.append(to_ms(mm))
vars_chart = [str(v) for v in V if str(v) != chart]
fname = 'hessrank_%s.ms'%chart
with open(fname,'w') as f:
    f.write(','.join(vars_chart)+'\n'); f.write('%d\n'%p)
    f.write(',\n'.join(eqs_chart)+'\n')
print('wrote', fname, len(eqs_chart), 'equations', flush=True)
t1=time.time()
nthreads = int(os.environ.get('NTHREADS','16'))
r = subprocess.run(['msolve','-v','1','-t',str(nthreads),'-f',fname,'-o','hessrank_%s.out'%chart], capture_output=True, text=True, timeout=7000)
print(r.stdout[-3000:], r.stderr[-3000:], flush=True)
out = open('hessrank_%s.out'%chart).read().strip()
print('OUTPUT:', out[:2000], flush=True)
res = {"exp":"J1-hessrank2-locus","family":"points v with rank Hess F(v)<=2","chart":chart,"d":None,"prime":p,
       "result": ("empty (GB={1})" if out.startswith('[1') or out.strip()=='[]' or 'no solution' in r.stdout.lower() else "see hessrank_%s.out"%chart),
       "raw_head": out[:200], "seconds": round(time.time()-t1,1), "job":"construct-dp-hessrank-%s"%chart,
       "note":"rank<=2 Hessian points: necessary for §2.7 exchange covers (rank 1, v not in X) and for cone tangent-hyperplane sections (v in X)"}
with open('wave2-construct.log.jsonl','a') as f: f.write(json.dumps(res)+'\n')
print(res)
