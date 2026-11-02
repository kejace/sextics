import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from w2lib import _col
tag, par = sys.argv[1], sys.argv[2]
pD = None if par == 'mixed' else int(par)
rr = [json.loads(l) for l in open(f'{WD}/{tag}.jsonl')]
rs = sorted([r for r in rr if 'seed' in r and 'error' not in r and 'a' in r], key=lambda r: r['r_part_rel'])
best = rs[0]; print('best stored:', best['seed'], best['r_part_rel'], 'status', best['status'])
fam = JacobianFamily(int(sys.argv[3]), pD, want_C=True, verbose=False)
a = np.array(best['a']); u0 = fam.P.T @ a; u0 /= np.linalg.norm(u0)
traj = []
def fun(u):
    r = fam.r_residual_u(u); traj.append(float(np.linalg.norm(r))); return r
t0 = time.time()
res = least_squares(fun, u0, jac=fam.r_jac_u, method='lm', max_nfev=int(sys.argv[4]) if len(sys.argv) > 4 else 600, xtol=1e-15, ftol=1e-15, gtol=1e-15)
a1, _ = fam.a_of_u(res.x)
D = _col(fam.Dfull @ a1); Om = fam.pat.omega(D); y0 = _col(fam.Nfull @ a1); b = Om @ y0
rp = float(np.linalg.norm(fam.C.T @ b)); bn = float(np.linalg.norm(b))
print('continuation: nfev', res.nfev, 'status', res.status, 'r_part_rel', rp / bn, 'time', round(time.time() - t0, 1))
print('trajectory (||r||):', [round(traj[i] / bn, 6) for i in [0, 1, 2, 5, 10, 20, 50, 100, 200, len(traj) - 1] if i < len(traj)])
# gradient norm at the end (stationarity)
J = fam.r_jac_u(res.x); g = J.T @ fam.r_residual_u(res.x)
print('|grad| / |r| =', float(np.linalg.norm(g) / max(np.linalg.norm(fam.r_residual_u(res.x)), 1e-300)))
jlog(f'{WD}/wave2-search.log.jsonl', {'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'wave': 2, 'family': 'stage-1 continuation from best point', 'tag': tag, 'd': int(sys.argv[3]), 'par': par,
     'start_r_part_rel': best['r_part_rel'], 'end_r_part_rel': rp / bn, 'nfev': int(res.nfev), 'status': int(res.status)})
