"""Exact replay (over Q) of fqlx's weighted-K3 Gauss-Manin obstruction at (t1,t2) = (2,3) and at a second point (t1,t2) = (5,-7).
F_sub(X0,Y0,U,V;t1,t2) weighted (1,1,2,2) degree 6; R = Q[X0,Y0,U,V]/J(F_sub); checks: weighted Hilbert series of R,
dim S_12 = 140, dim J_12 = 139, socle functional lambda (lambda(U^6)=1), and lambda(A dF_sub/dt_i) != 0 for A = X0^3 Y0^3.
All linear algebra exact (flint fmpq_mat)."""
import sys, json, time, itertools
from fractions import Fraction
import sympy as sp
from flint import fmpq_mat, fmpq
X0, Y0, U, V, t1, t2 = sp.symbols('X0 Y0 U V t1 t2')
W = {X0: 1, Y0: 1, U: 2, V: 2}
Fsub = (X0**6 + 5*X0**4*Y0**2 + 5*X0**2*Y0**4 + Y0**6 + (1 + t1**6)*U**3 + (1 + t2**6)*V**3
        + 5*(X0**2 + Y0**2)*(t1**2*U**2 + t1*t2*U*V + t2**2*V**2) + sp.Rational(5, 7)*t1*(1 + 2*t2**2)*X0*Y0*U*V)
# consistency with the instance: substitute U = X1^2, V = X2^2, t1 = Y1/X1, t2 = Y2/X2 and compare with F
X1, X2, Y1, Y2 = sp.symbols('X1 X2 Y1 Y2')
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
F, A, q, rec = core.load_instance('quadric-contact-residual')
Fsym = sum(sp.Rational(c.numerator, c.denominator) * X0**e[0]*X1**e[1]*X2**e[2]*Y0**e[3]*Y1**e[4]*Y2**e[5] for e, c in F.items())
chk = sp.simplify(sp.expand(Fsub.subs({U: X1**2, V: X2**2, t1: Y1/X1, t2: Y2/X2})) - Fsym)
print('F_sub consistent with F:', chk == 0)
def wmons(d):
    out = []
    for a in range(d + 1):
        for b in range(d + 1 - a):
            r = d - a - b
            if r % 2: continue
            for u in range(r // 2 + 1):
                out.append((a, b, u, r // 2 - u))
    return out
def tovec(poly, mons):
    P = sp.Poly(sp.expand(poly), X0, Y0, U, V)
    d = dict(P.terms())
    return [sp.Rational(d.get(m, 0)) for m in mons]
results = {}
for pt in [(2, 3), (5, -7)]:
    Fp = Fsub.subs({t1: pt[0], t2: pt[1]})
    dF = [sp.diff(Fp, v) for v in (X0, Y0, U, V)]
    hilb = []
    for d in range(0, 14):
        mons = wmons(d); idx = {m: i for i, m in enumerate(mons)}
        rows = []
        for v, g in zip((X0, Y0, U, V), dF):
            dg = d - (6 - W[v])
            if dg < 0: continue
            for m in wmons(dg):
                rows.append(tovec(g * X0**m[0]*Y0**m[1]*U**m[2]*V**m[3], mons))
        if rows:
            M = fmpq_mat(len(rows), len(mons), [fmpq(x.p, x.q) for r in rows for x in r])
            rk = M.rank()
        else:
            rk = 0
        hilb.append(len(mons) - rk)
        if d == 12:
            # socle functional: left nullspace of J_12 (rows) ... we need lambda with lambda(J_12)=0, lambda(U^6)=1
            # lambda is a vector in Q^140 with M lambda = 0 (M rows = generators); nullspace of M
            Mt = M.transpose()  # 140 x nrows
            # solve: find vector l with M l = 0  -> nullspace of M (as a map Q^140 -> Q^nrows)
            R, rk2 = M.rref()
            piv = []
            for i in range(rk2):
                j = 0
                while R[i, j] == 0: j += 1
                piv.append(j)
            free = [j for j in range(len(mons)) if j not in set(piv)]
            assert len(free) == 1, free
            f = free[0]
            lam = [fmpq(0)] * len(mons); lam[f] = fmpq(1)
            for i, j in enumerate(piv): lam[j] = -R[i, f]
            assert all(sum(M[r, c] * lam[c] for c in range(len(mons))) == 0 for r in range(0, M.nrows(), 37))
            u6 = lam[idx[(0, 0, 6, 0)]]
            assert u6 != 0
            lam = [x / u6 for x in lam]
            Apoly = X0**3 * Y0**3
            vals = {}
            for nm, dt in (('t1', sp.diff(Fsub, t1)), ('t2', sp.diff(Fsub, t2))):
                g = sp.expand(Apoly * dt.subs({t1: pt[0], t2: pt[1]}))
                vec = tovec(g, mons)
                vals[nm] = str(sum(fmpq(x.p, x.q) * l for x, l in zip(vec, lam)))
            results[str(pt)] = {'dimS12': len(mons), 'dimJ12': rk, 'lambda_A_dF_dt1': vals['t1'], 'lambda_A_dF_dt2': vals['t2']}
    results[str(pt)]['hilbert'] = hilb
    print(pt, results[str(pt)])
rec = {'exp': 'k3-gauss-manin-replay', 'prime': None, 'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'results': results,
       'agrees_with_fqlx': results['(2, 3)']['hilbert'] == [1, 2, 5, 8, 12, 14, 16, 14, 12, 8, 5, 2, 1, 0] and results['(2, 3)']['dimJ12'] == 139
                           and results['(2, 3)']['lambda_A_dF_dt1'] != '0' and results['(2, 3)']['lambda_A_dF_dt2'] != '0',
       'note': 'exact over Q (flint fmpq_mat); socle functional normalized lambda(U^6)=1; nonzero pairing => A dF_sub/dt_i not in J_12(F_sub) => residual class not GM-flat => no (t1,t2)-passive certificate at any degree (fqlx sec. 4)'}
open('/tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl', 'a').write(json.dumps(rec) + '\n')
print(json.dumps(rec))
