import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
from kcc import *
# validate the (unrestricted) WSystem against accepted witnesses, and the column construction against core.defect
for iid in ['cone-ladder-linear', 'two-reflection-cubic-surface-d2', 'f3-block-x0-quartic-strict', 'hesse-radial-hesse-cofactor-d4-rational']:
    F, A, q, rec = core.load_instance(iid)
    d, D, N, Q = core.load_witness(iid)
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    t = time.time()
    S = WSystem(F, A, q, D, P0, weights=False)
    r = S.solve(compress=False)
    # check: plug the known witness into the system columns: sum col_j x_j should equal rhs
    p = P0
    x = {}
    for (kind, i, m) in S.cols: pass
    acc = {}
    for j, (kind, i, m) in enumerate(S.cols):
        v = (N[i].get(m, 0) if kind == 'N' else Q.get(m, 0))
        v = core.reduce_mod({m: v}, p).get(m, 0) if v else 0
        if v:
            for e, c in S.colpolys[j].items():
                acc[e] = (acc.get(e, 0) + c * v) % p
    acc = {e: c for e, c in acc.items() if c}
    diff = core.sub(acc, S.rhs, p)
    print(f'{iid}: witness ok={ok}  d={d}  system rows={S.nrows} cols={S.ncols} consistent={r["consistent"]} rank={r["rank"]} nullity={r["nullity"]}  witness-satisfies-system={not diff}  [{time.time()-t:.1f}s]')
# also: a random-N defect check of the column construction on OUR F with a sigma-eigen D and weight restriction:
F, A, q = load()
import random
random.seed(1)
D = core.parse('X0**5*X1 + X1**5*X2 + X0*X2**5 + 3*X0*X1**3*Y0*Y1')  # weight 0, Y-even, d=6
S = WSystem(F, A, q, D, P0, weights=True)
# random N in the restricted monomial set, Q = 0; compare sum columns with core.defect - A D^2 (i.e. E - A D^2)
N = [{} for _ in range(NV)]
for (kind, i, m) in S.cols:
    if kind == 'N' and random.random() < 0.3:
        N[i][m] = random.randint(1, 50)
E = core.defect(F, A, D, q, N, P0)
AD2 = core.mul(A, core.mul(D, D, P0), P0)
lhs = core.sub(E, AD2, P0)   # = -q D N(F) + F(D div N - N(D))
acc = {}
for j, (kind, i, m) in enumerate(S.cols):
    v = N[i].get(m, 0) if kind == 'N' else 0
    if v:
        for e, c in S.colpolys[j].items():
            acc[e] = (acc.get(e, 0) + c * v) % P0
acc = {e: c for e, c in acc.items() if c}
print('column construction matches core.defect on our F (weight-restricted):', not core.sub(acc, lhs, P0))
print('restricted system size d=6:', S.nrows, S.ncols, 'unrestricted:', WSystem(F, A, q, D, P0, weights=False).ncols)
