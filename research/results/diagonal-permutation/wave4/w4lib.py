"""Wave-4 diagonal-permutation / second-reflection: shared helpers (msolve over F_p, logging)."""
import sys, json, time, subprocess, os, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import sympy as sp
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
Z6 = [X0,X1,X2,Y0,Y1,Y2]
p_,q_,a_,b_,t_ = sp.symbols('p q a b t'); W5 = [p_,q_,a_,b_,t_]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.expand(sp.sympify(core.to_str(F)))
LOG = '/tmp/hodge/results/diagonal-permutation/wave4/log.jsonl'
DIR = '/tmp/hodge/results/diagonal-permutation/wave4/'
P1 = 1073741827; P2 = 2147483647

def fibre(mu):
    """F restricted to the (01)-pencil member {X0-X1 = mu (Y0-Y1)} in coordinates (p,q,a,b,t):
    X0=(p+mu t)/2, X1=(p-mu t)/2, X2=a, Y0=(q+t)/2, Y1=(q-t)/2, Y2=b.  Even in t (reflection s)."""
    mu = sp.Rational(mu)
    Fm = sp.expand(Fs.subs({X0:(p_+mu*t_)/2, X1:(p_-mu*t_)/2, X2:a_, Y0:(q_+t_)/2, Y1:(q_-t_)/2, Y2:b_}, simultaneous=True))
    assert sp.expand(Fm.subs(t_, -t_) - Fm) == 0
    return Fm

def log(d):
    d = dict(d); d.setdefault('utc', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    d.setdefault('lane','diagonal-permutation'); d.setdefault('route','second-reflection')
    with open(LOG,'a') as f: f.write(json.dumps(d, default=str)+'\n')

def to_ms(poly, unknowns, p):
    poly = sp.Poly(sp.expand(poly), *unknowns)
    den = 1
    for c in poly.coeffs(): den = sp.ilcm(den, sp.Rational(c).q)
    poly = sp.Poly(sp.expand(poly.as_expr()*den), *unknowns)
    terms = []
    for mon, c in poly.terms():
        c = int(c) % p
        if c == 0: continue
        mstr = '*'.join('%s^%d'%(str(unknowns[i]),e) for i,e in enumerate(mon) if e>0)
        terms.append(('%d*%s'%(c, mstr)) if mstr else '%d'%c)
    return '+'.join(terms) if terms else '0'

def msolve(fname, unknowns, eqs, p, nthreads=2, timeout=3000):
    eqs = [sp.expand(e) for e in eqs]; eqs = [e for e in eqs if e != 0]
    with open(fname,'w') as f:
        f.write(','.join(str(x) for x in unknowns)+'\n'); f.write('%d\n'%p)
        f.write(',\n'.join(to_ms(e, unknowns, p) for e in eqs)+'\n')
    t1 = time.time()
    r = subprocess.run(['msolve','-v','0','-t',str(nthreads),'-f',fname,'-o',fname.replace('.ms','.out')],
                       capture_output=True, text=True, timeout=timeout)
    out = open(fname.replace('.ms','.out')).read().strip()
    return out, round(time.time()-t1,1), len(eqs)

def classify(out):
    s = out.strip()
    if s.startswith('[-1]'): return 'empty'
    if s.startswith('[0'): return 'finite'
    if s.startswith('[1'): return 'positive-dim'
    return 'unknown'

def parse_param(out):
    """Parse msolve zero-dim output -> (p, nvars, deg, vars, linform, elim, denom, params)."""
    import ast
    s = out.strip().rstrip(':')
    data = ast.literal_eval(s)
    assert data[0] == 0
    p, nvars, deg, vars_, lin, par = data[1][0], data[1][1], data[1][2], data[1][3], data[1][4], data[1][5]
    elim = par[1][0][1]; denom = par[1][1][1]; params = [pp[0][1] for pp in par[1][2]]
    return p, nvars, deg, vars_, lin, elim, denom, params

def solutions_mod_p(out):
    """Return list of F_p-rational solutions (as dict var->int) and the elim-poly degree / number of F_p roots."""
    import flint
    p, nvars, deg, vars_, lin, elim, denom, params = parse_param(out)
    R = flint.nmod_poly([c % p for c in elim], p)
    roots = [int(r) for r, m in R.roots()] if hasattr(R, 'roots') else [int(r) for r in R.factor()[1] and []]
    sols = []
    for T in roots:
        den = flint.nmod_poly([c % p for c in denom], p)(T)
        if int(den) == 0: sols.append(None); continue
        vals = [(-int(flint.nmod_poly([c % p for c in pp], p)(T))) * pow(int(den), -1, p) % p for pp in params]
        vals.append(T)   # last variable = linear form value (msolve convention when lin form is a new variable)
        sols.append(dict(zip(vars_, vals)))
    return {'deg_ideal': deg, 'elim_degree': len(elim)-1, 'n_roots_Fp': len(roots), 'vars': vars_, 'lin': lin, 'sols': sols}
