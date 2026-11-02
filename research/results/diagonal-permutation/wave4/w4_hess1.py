"""Rank-1 Hessian locus (2x2 minors) of F on P^5 or of the fibre F_mu on P^4, over F_p, all affine charts.
For chart z_i = 1 we run (a) the plain minor system and (b) for each other variable z_j the system
{minors, u*z_j - 1}; (b) all empty  <=>  the only F_p-bar point in the chart is the coordinate point e_i.
Usage: python3 w4_hess1.py P5|fibre <mu> <prime> [nthreads]"""
import sys, itertools, json, time
import sympy as sp
import w4lib as L
target = sys.argv[1]; mu = sys.argv[2]; p = int(sys.argv[3]); nth = int(sys.argv[4]) if len(sys.argv)>4 else 2
if target == 'P5':
    G = L.Fs; V = L.Z6
else:
    G = L.fibre(mu); V = L.W5
n = len(V)
H = sp.Matrix(n, n, lambda i,j: sp.expand(sp.diff(G, V[i], V[j])))
minors = []
seen = set()
for rows in itertools.combinations(range(n),2):
    for cols in itertools.combinations(range(n),2):
        if rows <= cols:
            m = sp.expand(H.extract(list(rows), list(cols)).det())
            if m != 0: minors.append(m)
print(target, 'mu', mu, 'p', p, 'minors', len(minors), flush=True)
u = sp.Symbol('u')
summary = {}
for i in range(n):
    sub = {V[i]: 1}
    eqs_i = [m.subs(sub) for m in minors]
    others = [V[j] for j in range(n) if j != i]
    tag = '%s_%s_p%d_chart%s' % (target, mu.replace('/','over'), p, V[i])
    out, sec, ne = L.msolve(L.DIR+'hess1_%s_plain.ms'%tag, others, eqs_i, p, nth)
    cls = L.classify(out)
    info = {'plain': cls, 'seconds': sec}
    if cls == 'finite':
        try:
            sol = L.solutions_mod_p(out); info['plain_sols'] = sol
        except Exception as e: info['plain_parse_error'] = str(e); info['raw'] = out[:300]
    res_j = {}
    for zj in others:
        out2, sec2, _ = L.msolve(L.DIR+'hess1_%s_u%s.ms'%(tag, zj), others+[u], eqs_i+[u*zj-1], p, nth)
        res_j[str(zj)] = (L.classify(out2), sec2)
    info['u_tests'] = res_j
    only_coord = all(v[0]=='empty' for v in res_j.values())
    info['only_coordinate_point'] = only_coord
    summary[str(V[i])] = info
    print(V[i], json.dumps(info, default=str)[:600], flush=True)
    L.log({'what':'rank-1 Hessian locus (2x2 minors), chart %s=1'%V[i], 'target':target, 'mu':mu, 'prime':p,
           'n_minors':len(minors), 'result': info, 'file': 'hess1_%s_*.ms/.out'%tag})
def chart_ok(v):
    if not v['only_coordinate_point']: return False
    if v['plain'] == 'empty': return True
    if v['plain'] == 'finite' and 'plain_sols' in v:
        return all(all(val == 0 for val in s.values()) for s in v['plain_sols']['sols'] if s is not None)
    return False
allok = all(chart_ok(v) for v in summary.values())
print('ALL CHARTS: rank-1 locus = coordinate points only:', allok, flush=True)
L.log({'what':'rank-1 Hessian locus summary', 'target':target, 'mu':mu, 'prime':p, 'result':'coordinate points only' if allok else 'SEE CHARTS', 'summary':summary})
