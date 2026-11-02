import sys, time, json, pickle
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from charsolve import *
import invsolve
F, A, q = load(); p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
QX, QY = e_k(X, 2), e_k(Y, 2); e1X, e1Y = e_k(X, 1), e_k(Y, 1); e3X, e3Y = e_k(X, 3), e_k(Y, 3)
e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e1Y, e3Y), 2))
D7, D7Y = mul(DX, e2sqX), mul(DY, e2sqY)
GsX = Group([[0,1,2],[3,4,5]], [-1, 1]); GsY = Group([[0,1,2],[3,4,5]], [1, -1]); GG = Group([[0,1,2],[3,4,5]], [1,1])
Treps = invsolve.orbit_reps(6, 1, invsolve.canon_codes)
def vec(P): return [P.get(r, 0) % p for r in Treps]
from flint import nmod_mat
def rank_of(vs):
    if not vs: return 0
    return nmod_mat(len(vs), len(Treps), [v[i] for v in vs for i in range(len(Treps))], p).rank()
Jp = [core.reduce_mod(P, p) for P in J_char_basis(F, 6, 1, GG, p)]
Jv = [vec(P) for P in Jp]; Av = vec(core.reduce_mod(A, p))
print('rank J^G_6 odd:', rank_of(Jv))
span = list(Jv); served = {}
for lab, D, G_ in [('D7', D7, GsX), ('D7Y', D7Y, GsY)]:
    t0 = time.time()
    S = CharSystem(F, A, q, D, G_, p=p, targets=True, verbose=True)
    tg = S.served_targets()
    vs = [vec(P) for P in tg]
    served[lab] = vs
    beyond = rank_of(Jv + vs) - rank_of(Jv)
    new = rank_of(span + vs) - rank_of(span); span += vs
    inA = rank_of(span + [Av]) == rank_of(span)
    msg = f'served^G total {len(vs)}, beyond J {beyond}, new {new}, cumulative beyond J {rank_of(span)-rank_of(Jv)}, A in span: {inA}'
    print(lab, msg, f'({time.time()-t0:.1f}s)', flush=True)
    log(lab, core.degree(D), 'served-char', True, None, None, None, time.time() - t0, msg + ' engine=charsolve')
pickle.dump(served, open(f'{RES}/served_d7.pkl', 'wb'))
