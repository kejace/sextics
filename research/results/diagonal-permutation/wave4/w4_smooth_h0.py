"""(a) smoothness of G = F (P^5) or F_mu (P^4) over F_p-bar: Jacobian ideal empty in every affine chart (msolve);
(b) h^0(T_X) = dim{M in gl_n : xi_M(G) in k*G} - 1, computed mod p and over Q by linear algebra.
Usage: python3 w4_smooth_h0.py P5|fibre <mu> <prime> [nthreads]"""
import sys, json, itertools
import sympy as sp
import w4lib as L
target = sys.argv[1]; mu = sys.argv[2]; p = int(sys.argv[3]); nth = int(sys.argv[4]) if len(sys.argv)>4 else 2
if target == 'P5': G = L.Fs; V = L.Z6
else: G = L.fibre(mu); V = L.W5
n = len(V)
grads = [sp.expand(sp.diff(G, v)) for v in V]
smooth = {}
for i in range(n):
    eqs = [g.subs({V[i]:1}) for g in grads]
    others = [V[j] for j in range(n) if j != i]
    tag = 'jac_%s_%s_p%d_chart%s' % (target, mu.replace('/','over'), p, V[i])
    out, sec, ne = L.msolve(L.DIR+tag+'.ms', others, eqs, p, nth)
    smooth[str(V[i])] = (L.classify(out), sec)
    print('chart', V[i], smooth[str(V[i])], flush=True)
is_smooth = all(v[0]=='empty' for v in smooth.values())
L.log({'what':'smoothness over F_p-bar via Jacobian ideal in all charts','target':target,'mu':mu,'prime':p,'result':{'smooth':is_smooth,'charts':smooth}})
# (b) h^0(T)
monos = sorted(sp.itermonomials(V, 6, 6), key=str) if False else None
cols = {}
vecs = []
labels = []
for i in range(n):
    for j in range(n):
        xi = sp.expand(V[j]*grads[i])          # vector field z_j d/dz_i applied to G
        vecs.append(sp.Poly(xi, *V)); labels.append('E%d%d'%(i,j))
vecs.append(sp.Poly(G, *V)); labels.append('G')
allmon = set()
for P in vecs: allmon.update(P.monoms())
allmon = sorted(allmon)
idx = {mm:k for k,mm in enumerate(allmon)}
rows = []
for P in vecs:
    row = [0]*len(allmon)
    for mm, c in P.terms(): row[idx[mm]] = c
    rows.append(row)
Mq = sp.Matrix(rows)
den = 1
for c in Mq: den = sp.ilcm(den, sp.Rational(c).q)
Mz = (Mq*den).applyfunc(lambda c: sp.Integer(c))
import flint
Mf = flint.fmpz_mat([[int(c) for c in Mz.row(r)] for r in range(Mz.rows)])
rankQ = Mf.rank()
Mp = flint.nmod_mat([[int(c) % p for c in Mz.row(r)] for r in range(Mz.rows)], p)
rankp = Mp.rank()
kerQ = n*n + 1 - rankQ; kerp = n*n + 1 - rankp
res = {'n':n, 'rows':n*n+1, 'cols':len(allmon), 'rank_Q':rankQ, 'rank_Fp':rankp, 'kernel_dim_Q':kerQ, 'kernel_dim_Fp':kerp,
       'h0_T_Q': kerQ-1, 'h0_T_Fp': kerp-1}
print(json.dumps(res), flush=True)
L.log({'what':'h^0(T_X) = dim{M: xi_M G in k G} - 1 by linear algebra (Q and F_p)','target':target,'mu':mu,'prime':p,'result':res})
