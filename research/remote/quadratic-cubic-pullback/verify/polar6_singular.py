"""Export the degree-6 Hessian-relaxation quadrics (phi(w) in J_6, w linear, pure parity) to Singular mod p and test
whether the projective zero set is empty (affine dim 0)."""
import sys, subprocess
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
from flint import nmod_mat
p = 32003
F, A, q = V.load_instance_F_A()
Fp = V.modp(F, p)
Fi = [V.pderiv(Fp, i, p) for i in range(6)]
H = [[V.pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
mons6e = V.monomials(6, 0); idx = {e: j for j, e in enumerate(mons6e)}
J6e = [V.shift(Fi[i], m, 1, p) for i in range(6) for m in V.monomials(1, (0 - (0 if i < 3 else 1)) % 2)]
MJ = V.to_matrix_rows(J6e, mons6e, p)
Jr, rk = MJ.rref(); assert rk == 18
piv = []
for r in range(rk):
    j = 0
    while int(Jr[r, j]) == 0: j += 1
    piv.append(j)
nonpiv = [j for j in range(len(mons6e)) if j not in set(piv)]
def nf(P):
    v = [0]*len(mons6e)
    for e, c in P.items(): v[idx[e]] = c
    for r, j in enumerate(piv):
        c = v[j]
        if c:
            for k in range(len(mons6e)):
                x = int(Jr[r, k])
                if x: v[k] = (v[k] - c*x) % p
    return [v[j] for j in nonpiv]
for mode, pD in (('even', 0), ('odd', 1)):
    basis = [(i, m) for i in range(6) for m in V.monomials(1, (pD - (0 if i < 3 else 1)) % 2)]
    if pD == 0: basis = [b for b in basis if b != (0, (1, 0, 0, 0, 0, 0))]   # gauge fix w -> w + l Z
    n = len(basis)
    eqs = {}
    for a in range(n):
        for b in range(a, n):
            (i, mi), (j, mj) = basis[a], basis[b]
            P_ = V.shift(H[i][j], tuple(x + y for x, y in zip(mi, mj)), 1, p)
            v = nf(P_)
            for r, c in enumerate(v):
                if c: eqs.setdefault(r, []).append((c, a, b))
    lines = []
    for r, terms in eqs.items():
        lines.append('+'.join(f'{c}*c({a})*c({b})' for c, a, b in terms))
    sing = f'ring r = {p}, (c(0..{n-1})), dp;\nideal I = ' + ',\n'.join(lines) + ';\nideal G = std(I);\nprint("mode {mode}: nvars {n}, neqs {len(lines)}, dim = " + string(dim(G)) + ", vdim = " + string(vdim(G)));\nquit;\n'
    fn = f'/tmp/hodge/results/quadratic-cubic-pullback/verify/polar6_{mode}.sing'
    open(fn, 'w').write(sing)
    out = subprocess.run(['timeout', '600', 'Singular', '-q', fn], capture_output=True, text=True).stdout.strip()
    print(out)
    empty = ('dim = 0' in out)
    log(f'C28 no {mode} degree-6 polar denominator passes E2 (Hessian relaxation phi(w) in J_6; Singular mod {p}: affine dim of the quadric ideal in {n} gauge-fixed unknowns)',
        'CONFIRMED' if empty else 'UNVERIFIED', out + ' | dim 0 of a homogeneous ideal = projectively empty = only w=0 (gauge) = D multiple of F; mod-p dim 0 implies dim 0 over Q (ideal membership of powers lifts)')
