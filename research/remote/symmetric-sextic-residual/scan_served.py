"""Served invariant-target spaces (sgn⊗sgn residue sector, 17-dim beyond J) for structured
semi-invariant denominators; accumulate the span and test whether [A] enters."""
import sys, time, json, pickle, itertools
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from charsolve import *
import invsolve
from flint import nmod_mat
F, A, q = load(); p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
def dlt(vs): return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
e = {('X',k): e_k(X,k) for k in (1,2,3)}; e.update({('Y',k): e_k(Y,k) for k in (1,2,3)})
QX, QY, e1X, e1Y, e3X, e3Y = e[('X',2)], e[('Y',2)], e[('X',1)], e[('Y',1)], e[('X',3)], e[('Y',3)]
p2X, p2Y, p3X, p3Y, p4X, p4Y, p6X, p6Y = [p_k(V, k) for V in (X, Y) for k in (2,)] + [p_k(V, k) for V in (X, Y) for k in (3,)] + [p_k(V,4) for V in (X,Y)] + [p_k(V,6) for V in (X,Y)]
p2X, p2Y = p_k(X,2), p_k(Y,2); p3X, p3Y = p_k(X,3), p_k(Y,3); p4X, p4Y = p_k(X,4), p_k(Y,4); p6X, p6Y = p_k(X,6), p_k(Y,6)
Fi = [core.pderiv(F, i) for i in range(6)]
def radial(idx):
    out = {}
    for i in idx: out = add(out, mul(var(i), Fi[i]))
    return out
EXF, EYF = radial((0,1,2)), radial((3,4,5))
GG = Group([[0,1,2],[3,4,5]], [1,1]); GsX = Group([[0,1,2],[3,4,5]], [-1,1]); GsY = Group([[0,1,2],[3,4,5]], [1,-1]); Gss = Group([[0,1,2],[3,4,5]], [-1,-1])
groups = {'triv': GG, 'sgnX': GsX, 'sgnY': GsY, 'sgnsgn': Gss}
Treps = invsolve.orbit_reps(6, 1, invsolve.canon_codes)
def vec(P): return [P.get(r, 0) % p for r in Treps]
def rank_of(vs):
    if not vs: return 0
    return nmod_mat(len(vs), len(Treps), [v[i] for v in vs for i in range(len(Treps))], p).rank()
Jv = [vec(core.reduce_mod(P, p)) for P in J_char_basis(F, 6, 1, GG, p)]
Av = vec(core.reduce_mod(A, p)); rJ = rank_of(Jv)
span = list(Jv); results = []
def typ(D):
    for lab, G_ in groups.items():
        if G_.is_semi_invariant(D): return lab, G_
    return None, None
def run(lab, D, fam):
    global span
    t0 = time.time()
    tl, G_ = typ(D)
    if tl is None or core.parity_of(D) is None:
        print(f'{lab:30s} skipped (not semi-invariant or mixed parity)', flush=True); return
    d = core.degree(D)
    S = CharSystem(F, A, q, D, G_, p=p, targets=True)
    tg = S.served_targets()
    vs = [vec(P) for P in tg]
    beyond = rank_of(Jv + vs) - rJ
    new = rank_of(span + vs) - rank_of(span); span = span + vs
    cum = rank_of(span) - rJ
    inA = rank_of(span + [Av]) == rank_of(span)
    # also consistency for A itself
    r = S.solve()
    msg = f'type={tl} served^G={len(vs)} beyondJ={beyond} new={new} cum_beyondJ={cum} A_in_span={inA} | A-system consistent={r["consistent"]} rank={r["rank"]} null={r["nullity"]} rows={S.nrows} cols={S.ncols-S.nT}'
    log(lab, d, fam, None, r['consistent'], r['rank'], r['nullity'], time.time() - t0, msg + ' engine=charsolve')
    print(f'{lab:30s} d={d} {msg} ({time.time()-t0:.1f}s)', flush=True)
    results.append((lab, d, tl, len(vs), beyond, new, cum, inA))
    if inA: print('!!!!! A ENTERS THE SPAN', flush=True)
    if beyond > 0: pickle.dump((lab, D, tg), open(f'{RES}/served_hit_{lab.replace(" ","_").replace("/","_").replace("*","x")}.pkl', 'wb'))
cands = []
C = lambda lab, D, fam: cands.append((lab, D, fam))
# sign-type low degree
C('DX', DX, 'sgn3'); C('DY', DY, 'sgn3')
for l, L in [('e1X', e1X), ('e1Y', e1Y)]:
    C(f'DX*{l}', mul(DX, L), 'sgn4'); C(f'DY*{l}', mul(DY, L), 'sgn4')
quads = [('QX', QX), ('p2X', p2X), ('QY', QY), ('p2Y', p2Y), ('e1X^2', pw(e1X,2)), ('e1Xe1Y', mul(e1X,e1Y)), ('e1Y^2', pw(e1Y,2)), ('QX-QY', sub(QX,QY)), ('QX+2QY', add(QX, sc(QY,2))), ('p2X+p2Y', add(p2X,p2Y))]
for l, Qd in quads:
    C(f'DX*{l}', mul(DX, Qd), 'sgn5'); C(f'DY*{l}', mul(DY, Qd), 'sgn5')
cubs = [('e3X', e3X), ('e1XQX', mul(e1X,QX)), ('e1X^3', pw(e1X,3)), ('p3X', p3X), ('e3Y', e3Y), ('e1YQY', mul(e1Y,QY)), ('e1Y^3', pw(e1Y,3)), ('p3Y', p3Y), ('e1XQY', mul(e1X,QY)), ('e1YQX', mul(e1Y,QX)), ('e1Xe1Y^2', mul(e1X, pw(e1Y,2))), ('e1X^2e1Y', mul(pw(e1X,2), e1Y))]
for l, Cb in cubs:
    C(f'DX*{l}', mul(DX, Cb), 'sgn6'); C(f'DY*{l}', mul(DY, Cb), 'sgn6')
C('DX*DY', mul(DX, DY), 'sgnsgn6'); C('DX*DY*e1X', mul(mul(DX, DY), e1X), 'sgnsgn7'); C('DX*DY*e1Y', mul(mul(DX, DY), e1Y), 'sgnsgn7')
for l, Qd in quads: C(f'DX*DY*{l}', mul(mul(DX, DY), Qd), 'sgnsgn8')
# invariant degree 6 structured (never run in the crashed inv_scan)
for a, b in [(1, -1), (1, 2), (2, 1), (1, -2), (-2, 1), (1, 3), (3, 1), (1, -3), (-3, 1), (1, 5), (5, 1), (2, -1), (-1, 2), (1, 1)]:
    C(f'{a}EX(F)+{b}EY(F)', add(sc(EXF, Fraction(a)), sc(EYF, Fraction(b))), 'J6-radial-pencil')
C('p6X-p6Y', sub(p6X, p6Y), 'fermat-nbr'); C('p6X+p6Y', add(p6X, p6Y), 'fermat-nbr')
C('DX^2+DY^2', add(pw(DX,2), pw(DY,2)), 'disc6'); C('DX^2-DY^2', sub(pw(DX,2), pw(DY,2)), 'disc6'); C('DX^2', pw(DX,2), 'disc6'); C('DY^2', pw(DY,2), 'disc6')
C('e3X^2+e3Y^2', add(pw(e3X,2), pw(e3Y,2)), 'sym6'); C('p3X^2+p3Y^2', add(pw(p3X,2), pw(p3Y,2)), 'sym6'); C('p3X^2-p3Y^2', sub(pw(p3X,2), pw(p3Y,2)), 'sym6')
C('QX^3+QY^3', add(pw(QX,3), pw(QY,3)), 'sym6'); C('QX^3-QY^3', sub(pw(QX,3), pw(QY,3)), 'sym6'); C('QX*QY*(QX+2QY)=M', mul(mul(QX, QY), add(QX, sc(QY, 2))), 'coupling')
C('e3X*e3Y=A', mul(e3X, e3Y), 'sym6-odd'); C('p3X*p3Y', mul(p3X, p3Y), 'sym6-odd'); C('e3X*p3Y', mul(e3X, p3Y), 'sym6-odd'); C('p3X*e3Y', mul(p3X, e3Y), 'sym6-odd')
# quotient-curve family p6X + cubic(QX,QY), p6Y + cubic
vals = [0, 1, -1, 2, -2, Fraction(1,2), Fraction(-1,2)]
for al, be, ga, de in itertools.product(vals, repeat=4):
    if (al, be, ga, de) == (0,0,0,0): continue
    cub = add(add(sc(pw(QX,3), al), sc(mul(pw(QX,2), QY), be)), add(sc(mul(QX, pw(QY,2)), ga), sc(pw(QY,3), de)))
    C(f'p6X+({al},{be},{ga},{de})', add(p6X, cub), 'quot-curve6X')
if __name__ == '__main__':
    only = sys.argv[1:]
    for lab, D, fam in cands:
        if only and not any(o in fam for o in only): continue
        run(lab, D, fam)
    pickle.dump(results, open(f'{RES}/scan_served_results.pkl', 'wb'))
