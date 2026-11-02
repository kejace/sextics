"""Rank <= 2 Hessian locus of F on P^5 (3x3 minors), over F_p, all six affine charts (wave-2 item J1, never completed).
Reports dimension, ideal degree, elimination polynomial degree and its factorisation pattern over F_p, and F_p-rational points.
Usage: python3 w4_hess2_p5.py <prime> [nthreads]"""
import sys, itertools, json, time
import sympy as sp, flint
import w4lib as L
p = int(sys.argv[1]); nth = int(sys.argv[2]) if len(sys.argv)>2 else 2
V = L.Z6; n = 6
H = sp.Matrix(n, n, lambda i,j: sp.expand(sp.diff(L.Fs, V[i], V[j])))
minors = []
for rows in itertools.combinations(range(n),3):
    for cols in itertools.combinations(range(n),3):
        if rows <= cols:
            m = sp.expand(H.extract(list(rows), list(cols)).det())
            if m != 0: minors.append(m)
print('3x3 minors', len(minors), flush=True)
for i in range(n):
    eqs = [m.subs({V[i]:1}) for m in minors]
    others = [V[j] for j in range(n) if j != i]
    tag = 'hess2_P5_p%d_chart%s' % (p, V[i])
    out, sec, ne = L.msolve(L.DIR+tag+'.ms', others, eqs, p, nth, timeout=5000)
    cls = L.classify(out)
    info = {'class': cls, 'seconds': sec}
    if cls == 'finite':
        pp, nv, deg, vars_, lin, elim, denom, params = L.parse_param(out)
        R = flint.nmod_poly([c % p for c in elim], p)
        fac = R.factor()
        info.update({'deg_ideal': deg, 'elim_degree': R.degree(), 'elim_factor_degrees': sorted(int(f.degree()) for f, e in fac[1]),
                     'elim_factor_mults': [int(e) for f, e in fac[1]]})
        sol = L.solutions_mod_p(out); info['Fp_points'] = sol['sols']
    else:
        info['raw'] = out[:300]
    print(V[i], json.dumps(info, default=str)[:1500], flush=True)
    L.log({'what':'rank<=2 Hessian locus of F on P^5 (3x3 minors), chart %s=1'%V[i], 'prime':p, 'n_minors':len(minors), 'result':info, 'file':tag+'.ms/.out'})
