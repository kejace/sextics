"""Exact (over Q) check that the rank-1 Hessian locus of F in chart z_i = 1 is only the coordinate point e_i:
Singular std over Q of the 2x2 minors ideal, then verify x_j^k in I for every other variable (radical = maximal ideal).
Charts X0 and Y0 suffice by the S3 x <tau> symmetry (S3 permutes X-charts and Y-charts, tau preserves them).
Usage: python3 w4_exactQ_hess1.py P5|fibre <mu> <chart>"""
import sys, itertools, subprocess, time, json
import sympy as sp
import w4lib as L
target, mu, chart = sys.argv[1], sys.argv[2], sys.argv[3]
if target == 'P5': G = L.Fs; V = L.Z6
else: G = L.fibre(mu); V = L.W5
n = len(V); ci = [str(v) for v in V].index(chart)
H = sp.Matrix(n, n, lambda i,j: sp.expand(sp.diff(G, V[i], V[j])))
minors = []
for rows in itertools.combinations(range(n),2):
    for cols in itertools.combinations(range(n),2):
        if rows <= cols:
            m = sp.expand(H.extract(list(rows), list(cols)).det().subs({V[ci]:1}))
            if m != 0: minors.append(m)
others = [V[j] for j in range(n) if j != ci]
def sing(e):
    P = sp.Poly(e, *others); den = 1
    for c in P.coeffs(): den = sp.ilcm(den, sp.Rational(c).q)
    return str(sp.expand(e*den)).replace('**','^')
tag = 'exactQ_hess1_%s_%s_%s' % (target, mu.replace('/','over'), chart)
with open(tag+'.sing','w') as f:
    f.write('ring r = 0,(%s),dp;\n' % ','.join(str(v) for v in others))
    f.write('ideal I = %s;\n' % ',\n'.join(sing(m) for m in minors))
    f.write('option(redSB); ideal G = std(I);\n')
    f.write('print("DIM"); print(dim(G)); print("VDIM"); print(vdim(G));\n')
    for v in others:
        f.write('int k%s = 1; while (k%s < 60 and reduce(%s^k%s, G) != 0) { k%s = k%s + 1; } print("POWER %s"); print(k%s);\n' % ((v,)*8))
    f.write('quit;\n')
t0 = time.time()
r = subprocess.run(['Singular','-q','--no-warn', tag+'.sing'], capture_output=True, text=True, timeout=5400)
out = r.stdout
print(out[-1500:], r.stderr[-500:], 'seconds', round(time.time()-t0,1), flush=True)
L.log({'what':'EXACT over Q (Singular std): rank-1 Hessian locus in chart %s=1; dim/vdim and minimal k with z_j^k in I'%chart,
       'target':target,'mu':mu,'prime':None,'result':out.replace('\n',' ')[-800:],'seconds':round(time.time()-t0,1),'file':tag+'.sing'})
