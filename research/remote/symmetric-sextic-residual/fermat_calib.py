import sys, time
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from charsolve import *
import invsolve
from flint import nmod_mat
F, A, q = load(); p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
F0 = add(p_k(X, 6), p_k(Y, 6))
GG = Group([[0,1,2],[3,4,5]], [1,1])
Treps = invsolve.orbit_reps(6, 1, invsolve.canon_codes)
def vec(P): return [P.get(r, 0) % p for r in Treps]
def rank_of(vs):
    if not vs: return 0
    return nmod_mat(len(vs), len(Treps), [v[i] for v in vs for i in range(len(Treps))], p).rank()
def served_inv(Fx, D, lab):
    S = CharSystem(Fx, A, q, D, GG, p=p, targets=True)
    tg = S.served_targets()
    J = [vec(core.reduce_mod(P, p)) for P in J_char_basis(Fx, 6, 1, GG, p)]
    vs = [vec(P) for P in tg]
    Av = vec(core.reduce_mod(A, p))
    print(f'{lab:28s} d={core.degree(D)} served^G={len(vs)} beyondJ={rank_of(J+vs)-rank_of(J)} A served: {rank_of(J+vs+[Av])==rank_of(J+vs)}', flush=True)
    return tg
print('=== Fermat F0 calibration (invariant served targets) ===')
served_inv(F0, p_k(X, 6), 'p6X (join, control)')
for c in (0, 1, -1, 2, -3, 6):
    served_inv(F0, add(p_k(X, 3), sc(e_k(X, 3), c)), f'Hesse p3X+({c})e3X')
served_inv(F0, e_k(X, 3), 'e3X (3 lines)')
served_inv(F0, p_k(X, 4), 'p4X (Fermat quartic cone)')
served_inv(F0, add(p_k(X, 4), sc(mul(e_k(X,2), e_k(X,2)), 1)), 'p4X+QX^2')
served_inv(F0, p_k(X, 5), 'p5X')
served_inv(F0, mul(e_k(X,1), p_k(X, 4)), 'e1X*p4X')
served_inv(F0, add(p_k(X, 5), mul(e_k(X,2), p_k(X,3))), 'p5X+QX*p3X')
served_inv(F0, add(p_k(X, 6), sc(pw(e_k(X,2),3), 1)), 'p6X+QX^3')
served_inv(F0, add(p_k(X, 6), sc(mul(e_k(X,3), p_k(X,3)), 1)), 'p6X+e3X*p3X')
served_inv(F0, add(p_k(X, 6), sc(p_k(Y,6), 1)), 'p6X+p6Y (=F0!)')
served_inv(F0, add(p_k(X, 3), p_k(Y, 3)), 'p3X+p3Y')
served_inv(F0, mul(p_k(X, 3), p_k(Y, 3)), 'p3X*p3Y')
print('=== our F: does D7 kill all of R_6^odd? ===')
DX = mul(mul(sub(X[0], X[1]), sub(X[0], X[2])), sub(X[1], X[2]))
QX = e_k(X, 2); e2sqX = sub(pw(QX, 2), sc(mul(e_k(X,1), e_k(X,3)), 2)); D7 = mul(DX, e2sqX)
allodd = [{m: 1} for m in core.monomials(6, 1)]
t0 = time.time(); rk, _ = sector_kernel(F, q, D7, allodd); print('kernel of D7 on all of R_6^odd (incl J, dim J_6^odd=?):', rk, f'({time.time()-t0:.1f}s)')
rk1, _ = sector_kernel(F, q, core.parse('1'), allodd); print('kernel of 1 (= dim J_6^odd):', rk1)
