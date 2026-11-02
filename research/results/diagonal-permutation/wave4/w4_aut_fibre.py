"""Full linear automorphism group of the fibre F_mu over F_p-bar, GIVEN (checked separately) that the rank-1 Hessian
locus of F_mu is exactly {e_a, e_b}.  Then every g in Aut(F_mu) permutes {e_a,e_b}, hence preserves L = span(e_a,e_b)
and M = {a = b = 0} (intersection of the polar hyperplanes {a=0},{b=0}), so g = (g_L, g_M) is block diagonal with
g_L in Aut(a^6+b^6) = monomial 2x2 with 6th roots of unity.  Projectively 12 cases for g_L; for each we solve
F_mu o (g_L, g_M) = F_mu for g_M in GL_3 (9 unknowns) with msolve over F_p.
Usage: python3 w4_aut_fibre.py <mu> <prime> [nthreads]"""
import sys, json, time, itertools
import sympy as sp
import w4lib as L
mu = sys.argv[1]; p = int(sys.argv[2]); nth = int(sys.argv[3]) if len(sys.argv)>3 else 2
assert p % 6 == 1
# primitive 6th root of unity mod p
g0 = 2
while pow(g0, (p-1)//2, p) == 1 or pow(g0, (p-1)//3, p) == 1: g0 += 1
zeta = pow(g0, (p-1)//6, p)
assert pow(zeta, 6, p) == 1 and pow(zeta, 3, p) != 1 and pow(zeta, 2, p) != 1
Fm = L.fibre(mu)
pv, qv, av, bv, tv = L.W5
m = sp.symbols('m0:9'); Mmat = sp.Matrix(3,3,m)
def subs_for(case, j):
    zj = pow(zeta, j, p)
    newp, newq, newt = (Mmat*sp.Matrix([pv,qv,tv]))
    d = {pv:newp, qv:newq, tv:newt}
    if case == 'diag': d.update({av: av, bv: zj*bv})
    else:               d.update({av: bv, bv: zj*av})
    return d
results = []
total = 0
for case in ['diag','swap']:
    for j in range(6):
        d = subs_for(case, j)
        G = sp.expand(Fm.subs(d, simultaneous=True) - Fm)
        P = sp.Poly(G, *L.W5)
        eqs = [sp.expand(c) for c in P.coeffs()]
        tag = 'autfib_%s_p%d_%s%d' % (mu.replace('/','over'), p, case, j)
        out, sec, ne = L.msolve(L.DIR+tag+'.ms', list(m), eqs, p, nth, timeout=6000)
        cls = L.classify(out)
        info = {'case':case, 'j':j, 'class':cls, 'seconds':sec, 'n_eqs':ne}
        if cls == 'finite':
            sol = L.solutions_mod_p(out)
            info['deg_ideal'] = sol['deg_ideal']; info['elim_degree'] = sol['elim_degree']; info['n_roots_Fp'] = sol['n_roots_Fp']
            mats = []
            for s in sol['sols']:
                if s is None: mats.append(None); continue
                vals = [s[str(mi)] for mi in m]
                # verify by substitution mod p
                dd = dict(subs_for(case, j)); Mnum = sp.Matrix(3,3,[sp.Integer(v) for v in vals])
                newp, newq, newt = (Mnum*sp.Matrix([pv,qv,tv]))
                dd.update({pv:newp, qv:newq, tv:newt})
                chk = sp.Poly(sp.expand(Fm.subs(dd, simultaneous=True) - Fm), *L.W5)
                den = 1
                for c in chk.coeffs(): den = sp.ilcm(den, sp.Rational(c).q)
                ok = all(int(c*den) % p == 0 for c in chk.coeffs())
                mats.append({'g_M_rows': [[(v if v <= p//2 else v-p) for v in vals[3*r:3*r+3]] for r in range(3)], 'verified_mod_p': ok})
            info['solutions'] = mats
            total += sol['deg_ideal']
        elif cls != 'empty':
            info['raw'] = out[:400]
        print(json.dumps(info, default=str), flush=True)
        L.log({'what':'fibre automorphisms, block case g_L=%s zeta^%d, solve for g_M (9 unknowns)'%(case,j), 'mu':mu,'prime':p,
               'zeta':zeta, 'result':info, 'file':tag+'.ms/.out'})
        results.append(info)
print('TOTAL number of projective linear automorphisms of F_mu over F_p-bar (sum of ideal degrees):', total, flush=True)
L.log({'what':'fibre automorphism group summary (given rank-1 locus = {e_a,e_b})','mu':mu,'prime':p,'result':{'order':total,'cases':results}})
