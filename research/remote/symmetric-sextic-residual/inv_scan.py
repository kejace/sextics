"""Invariant full solves (own reduced solver, engine-independent) on invariant denominators."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
QX, QY = e_k(X, 2), e_k(Y, 2); e1X, e1Y = e_k(X, 1), e_k(Y, 1); e3X, e3Y = e_k(X, 3), e_k(Y, 3)
p3X, p3Y = p_k(X, 3), p_k(Y, 3); p6X, p6Y = p_k(X, 6), p_k(Y, 6)
def dlt(vs):
    return mul(mul(sub(vs[0], vs[1]), sub(vs[0], vs[2])), sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
def radial(F, idx):
    out = {}
    for i in idx:
        out = add(out, mul(var(i), core.pderiv(F, i)))
    return out
EXF, EYF = radial(F, (0, 1, 2)), radial(F, (3, 4, 5))
e2sqX = sub(pw(QX, 2), sc(mul(e1X, e3X), 2)); e2sqY = sub(pw(QY, 2), sc(mul(e1Y, e3Y), 2))
D7, D7Y = mul(DX, e2sqX), mul(DY, e2sqY)
cands = []
def C(lab, D, fam): cands.append((lab, D, fam))
C('EX(F)', EXF, 'J6-radial'); C('EY(F)', EYF, 'J6-radial')
for a, b in [(1, -1), (1, 2), (2, 1), (1, -2), (-2, 1), (1, 3), (3, 1), (1, -3), (-3, 1), (1, 5), (5, 1), (2, -1), (-1, 2)]:
    C(f'{a}EX(F)+{b}EY(F)', add(sc(EXF, Fraction(a)), sc(EYF, Fraction(b))), 'J6-radial-pencil')
C('p6X', p6X, 'fermat-nbr'); C('p6Y', p6Y, 'fermat-nbr'); C('p6X+p6Y', add(p6X, p6Y), 'fermat-nbr'); C('p6X-p6Y', sub(p6X, p6Y), 'fermat-nbr')
C('DX^2+DY^2', add(pw(DX, 2), pw(DY, 2)), 'disc6'); C('DX^2-DY^2', sub(pw(DX, 2), pw(DY, 2)), 'disc6')
C('e3X^2+e3Y^2', add(pw(e3X, 2), pw(e3Y, 2)), 'sym6'); C('p3X^2+p3Y^2', add(pw(p3X, 2), pw(p3Y, 2)), 'sym6'); C('p3X^2-p3Y^2', sub(pw(p3X, 2), pw(p3Y, 2)), 'sym6')
C('QX^3+QY^3', add(pw(QX, 3), pw(QY, 3)), 'sym6'); C('QX^3-QY^3', sub(pw(QX, 3), pw(QY, 3)), 'sym6')
C('e3X*e3Y=A', mul(e3X, e3Y), 'sym6-odd'); C('p3X*p3Y', mul(p3X, p3Y), 'sym6-odd'); C('DX*DY', mul(DX, DY), 'disc6-odd')
C('D7*D7Y', mul(D7, D7Y), 'e1-kernel-product-d14')
if __name__ == '__main__':
    only = sys.argv[1:] 
    for lab, D, fam in cands:
        if only and lab not in only: continue
        d = core.degree(D)
        t0 = time.time()
        if not is_invariant(D):
            print(lab, 'NOT invariant, skipping'); continue
        ok, info = e1_inv(F, A, q, D)
        if not ok:
            log(lab, d, fam, False, None, info['rank'], info['cols'] - info['rank'], time.time() - t0, f'inv-solver E1 only: {info} engine=inv-v1')
            print(f'{lab:22s} d={d} E1=False {info} ({time.time()-t0:.1f}s)', flush=True); continue
        st, data = solve_inv(F, A, q, D, verbose=False)
        cons = st != 'inconsistent'
        log(lab, d, fam, True, cons, data.get('rank'), data.get('nullity'), time.time() - t0,
            f'inv-solver full: status={st} rows={data.get("rows")} cols={data.get("cols")} engine=inv-v1')
        print(f'{lab:22s} d={d} E1=True status={st} rank={data.get("rank")} null={data.get("nullity")} rows={data.get("rows")} cols={data.get("cols")} ({time.time()-t0:.1f}s)', flush=True)
        if st == 'solved':
            print('!!!!! SOLVED', lab, flush=True)
            json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}, open(f'{RES}/witness-candidate.json', 'w'), indent=1)
        elif cons:
            print('!!!!! consistent mod p but lift failed', lab, flush=True)
