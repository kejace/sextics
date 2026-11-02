"""Positive control of the sector VarPro pipeline on solved instances with explicit witnesses (q=1, monomial A).
Characters are taken from the stabilizer subgroup of (A, D_witness) so that the witness is an eigenvector."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual')
import fplib, numpy as np
from hodge import core, symmetry
iid = sys.argv[1]
nrest = int(sys.argv[2]) if len(sys.argv) > 2 else 4
F, A, q, G = fplib.load(iid)
w = json.load(open(f'/home/kejace/devel/sextics/data/{iid}/witness.json'))
Dw = core.parse(w['D']); d = core.degree(Dw)
sub = symmetry.stabilizer(G, [A, Dw])
ch = fplib.Chars(G, elements=sub)
c = ch.char_of(Dw)
print(iid, 'q', q, 'group', G, 'stabilizer size', len(sub), 'A', core.to_str(A), 'cA', ch.char_of(A), 'witness d', d, 'char', c)
dm = fplib.sector_dims_modp(F, A, q, d, c, ch); print('modp dims', dm)
B, Dm, info = fplib.sector_basis_numeric(F, A, q, d, c, ch, verbose=True); print(info)
idx = {m: i for i, m in enumerate(Dm)}
v = fplib.poly_to_vec(Dw, idx, len(Dm))
proj = B.T @ (B @ v); print('witness in span: rel err', np.linalg.norm(v - proj) / np.linalg.norm(v))
fam = fplib.SectorVarPro(F, A, q, B, Dm, c, ch, verbose=True, tag=iid)
print('fd check', fam.fd_check(0), fam.fd_check(1))
a = B @ v; u = fam.P.T @ a
print('relres at witness', fam.relres(u))
r = fam.lm_from(0, init=u); print('LM from witness', {k: r[k] for k in ('relres', 'iters', 'status', 'time')})
hits = 0
for s in range(nrest):
    r = fam.minimize_from(s, maxiter=300); print('LBFGS seed', s, {k: r[k] for k in ('relres', 'iters', 'time')})
    r2 = fam.lm_from(s, init=r['u']); print('   LM polish', {k: r2[k] for k in ('relres', 'iters', 'status', 'time')})
    if r2['relres'] < 1e-8:
        print('   D ~', fam.D_float_str(r2['Dc'])[:200])
        recs = fplib.try_reconstruct_and_confirm(fam, r2['Dc'], G=None)
        for rr in recs:
            print('   recon', {k: rr[k] for k in rr if k not in ('N', 'Q')})
        if any(rr.get('status') == 'solved' for rr in recs):
            hits += 1
            print('   CONFIRMED WITNESS'); break
print('hits', hits)
