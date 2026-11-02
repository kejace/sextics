"""Degree-7 Y-even polar relaxation: w_i quadratic (par 0 for i<3, 1 for i>=3), gauge l in S_1^even (3 dims), phi(w) in J_8^even
(valid since K_8 = J_8 both parities => A. injective on R_8). Writes Singular and msolve inputs; runs Singular (timeout from argv)."""
import sys, subprocess, time
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
from flint import nmod_mat
p = 32003
tmo = sys.argv[1] if len(sys.argv) > 1 else '3600'
F, A, q = V.load_instance_F_A()
Fp = V.modp(F, p)
Fi = [V.pderiv(Fp, i, p) for i in range(6)]
H = [[V.pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
mons = V.monomials(8, 0); idx = {e: j for j, e in enumerate(mons)}
J8 = [V.shift(Fi[i], m, 1, p) for i in range(6) for m in V.monomials(3, (0 - (0 if i < 3 else 1)) % 2)]
MJ = V.to_matrix_rows(J8, mons, p)
Jr, rk = MJ.rref(); assert rk == 168, rk
piv = []
for r in range(rk):
    j = 0
    while int(Jr[r, j]) == 0: j += 1
    piv.append(j)
pivset = set(piv); nonpiv = [j for j in range(len(mons)) if j not in pivset]
Jrows = [[int(Jr[r, k]) for k in range(len(mons))] for r in range(rk)]
def nf(P):
    v = [0]*len(mons)
    for e, c in P.items(): v[idx[e]] = c
    for r, j in enumerate(piv):
        c = v[j]
        if c:
            row = Jrows[r]
            for k in range(len(mons)):
                if row[k]: v[k] = (v[k] - c*row[k]) % p
    return [v[j] for j in nonpiv]
basis = [(i, m) for i in range(6) for m in V.monomials(2, (0 - (0 if i < 3 else 1)) % 2)]
# gauge: w -> w + l Z, l in S_1^even = {X0, X1, X2}: fix by dropping (0, X0*X_k) for k = 0,1,2
drop = {(0, (2,0,0,0,0,0)), (0, (1,1,0,0,0,0)), (0, (1,0,1,0,0,0))}
basis = [b for b in basis if b not in drop]
n = len(basis)
t0 = time.time()
eqs = {}
for a in range(n):
    for b in range(a, n):
        (i, mi), (j, mj) = basis[a], basis[b]
        P_ = V.shift(H[i][j], tuple(x + y for x, y in zip(mi, mj)), 1, p)
        for r, c in enumerate(nf(P_)):
            if c: eqs.setdefault(r, []).append((c, a, b))
lines = ['+'.join(f'{c}*c({a})*c({b})' for c, a, b in terms) for terms in eqs.values()]
print(f'd=7 even: {n} unknowns, {len(lines)} quadrics, built in {time.time()-t0:.1f}s', flush=True)
sing = f'ring r = {p}, (c(0..{n-1})), dp;\nideal I = ' + ',\n'.join(lines) + ';\nint t0 = timer;\nideal G = std(I);\nprint("d7 even: nvars {n}, neqs {len(lines)}, dim = " + string(dim(G)) + ", vdim = " + string(vdim(G)) + ", time " + string(timer - t0));\nquit;\n'
open('/tmp/hodge/results/quadratic-cubic-pullback/verify/polar7_even.sing', 'w').write(sing)
with open('/tmp/hodge/results/quadratic-cubic-pullback/verify/polar7_even.ms', 'w') as f:
    f.write(','.join(f'c{k}' for k in range(n)) + '\n' + f'{p}\n')
    f.write(',\n'.join(l.replace('c(', 'c').replace(')', '') for l in lines) + '\n')
out = subprocess.run(['timeout', tmo, 'Singular', '-q', '/tmp/hodge/results/quadratic-cubic-pullback/verify/polar7_even.sing'], capture_output=True, text=True).stdout.strip()
print(out, flush=True)
log('C29 no Y-even degree-7 polar denominator passes E2 (Hessian relaxation phi(w) in J_8^even, 60 gauge-fixed unknowns; Singular mod 32003 affine dim)',
    'CONFIRMED' if 'dim = 0' in out else 'UNVERIFIED', out or f'Singular timed out after {tmo}s')
