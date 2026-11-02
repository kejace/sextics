"""Low-degree family scan (degrees 1-4 full dense solve, degrees 5-6 compressed)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
F, A, q = load()
QX, QY = e_k(X, 2), e_k(Y, 2)
e1X, e1Y = e_k(X, 1), e_k(Y, 1)
e3X, e3Y = e_k(X, 3), e_k(Y, 3)
p2X, p2Y = p_k(X, 2), p_k(Y, 2)
p3X, p3Y = p_k(X, 3), p_k(Y, 3)
p4X, p4Y = p_k(X, 4), p_k(Y, 4)
p6X, p6Y = p_k(X, 6), p_k(Y, 6)
def dlt(vs):
    return core.mul(core.mul(core.sub(vs[0], vs[1]), core.sub(vs[0], vs[2])), core.sub(vs[1], vs[2]))
DX, DY = dlt(X), dlt(Y)
M = core.mul(core.mul(QX, QY), core.add(QX, core.scale(QY, 2)))   # mixed part /(5/7)
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
fams = []
def F_(label, D, fam):
    fams.append((label, D, fam))
# degree 1
F_('X0', X[0], 'linear'); F_('e1X', e1X, 'linear'); F_('X0-X1', sub(X[0], X[1]), 'linear'); F_('X0+X1', add(X[0], X[1]), 'linear')
F_('X0+X1-2X2', sub(add(X[0], X[1]), sc(X[2], 2)), 'linear')
F_('Y0', Y[0], 'linear'); F_('e1Y', e1Y, 'linear'); F_('Y0-Y1', sub(Y[0], Y[1]), 'linear')
# degree 2
F_('QX', QX, 'sym2'); F_('QY', QY, 'sym2'); F_('p2X', p2X, 'sym2'); F_('p2Y', p2Y, 'sym2')
for c in (1, -1, 2, -2, Fraction(1, 2), Fraction(-1, 2), 3, -3):
    F_(f'QX+({c})QY', add(QX, sc(QY, Fraction(c))), 'sym2-pencil')
for c in (1, -1):
    F_(f'p2X+({c})p2Y', add(p2X, sc(p2Y, Fraction(c))), 'sym2-pencil')
F_('e1X^2-e1Y^2', sub(pw(e1X, 2), pw(e1Y, 2)), 'sym2')
F_('X0^2-Y0^2', sub(pw(X[0], 2), pw(Y[0], 2)), 'quad'); F_('X0^2+Y0^2', add(pw(X[0], 2), pw(Y[0], 2)), 'quad')
F_('X0X1', mul(X[0], X[1]), 'quad'); F_('X0X1-Y0Y1', sub(mul(X[0], X[1]), mul(Y[0], Y[1])), 'quad')
F_('X0Y0', mul(X[0], Y[0]), 'quad-odd'); F_('e1Xe1Y', mul(e1X, e1Y), 'quad-odd'); F_('(X0-X1)(Y0-Y1)', mul(sub(X[0], X[1]), sub(Y[0], Y[1])), 'quad-odd')
F_('p2X-QX', sub(p2X, QX), 'sym2'); F_('p2Y-QY', sub(p2Y, QY), 'sym2')
# degree 3
for c in (0, 1, -1, 2, -2, 3, -3, 6, -6, Fraction(3, 2), Fraction(-3, 2)):
    F_(f'p3Y+({c})e3Y', add(p3Y, sc(e3Y, Fraction(c))), 'Ycubic-Hesse')
    F_(f'p3X+({c})e3X', add(p3X, sc(e3X, Fraction(c))), 'Xcubic-Hesse')
F_('e3Y', e3Y, 'Ycubic'); F_('e3X', e3X, 'Xcubic'); F_('e1Y*QY', mul(e1Y, QY), 'Ycubic'); F_('e1X*QX', mul(e1X, QX), 'Xcubic')
F_('DX', DX, 'disc'); F_('DY', DY, 'disc')
F_('e1X*QY', mul(e1X, QY), 'mixed3'); F_('e1X*p2Y', mul(e1X, p2Y), 'mixed3'); F_('e1Y*QX', mul(e1Y, QX), 'mixed3'); F_('e1Y*p2X', mul(e1Y, p2X), 'mixed3')
F_('X0*QY', mul(X[0], QY), 'mixed3'); F_('X0*p2Y', mul(X[0], p2Y), 'mixed3'); F_('X0Y0Y1', mul(X[0], mul(Y[0], Y[1])), 'mixed3')
F_('X0^3+X1^3+Y0^2*X2', add(add(pw(X[0],3), pw(X[1],3)), mul(pw(Y[0],2), X[2])), 'misc3')
# degree 4
F_('p4X', p4X, 'sym4'); F_('p4Y', p4Y, 'sym4'); F_('p4X+p4Y', add(p4X, p4Y), 'sym4'); F_('p4X-p4Y', sub(p4X, p4Y), 'sym4')
F_('QX*QY', mul(QX, QY), 'sym4'); F_('p2X*p2Y', mul(p2X, p2Y), 'sym4'); F_('QX*p2Y', mul(QX, p2Y), 'sym4'); F_('p2X*QY', mul(p2X, QY), 'sym4')
F_('QX^2-QY^2', sub(pw(QX, 2), pw(QY, 2)), 'sym4'); F_('p2X^2-p2Y^2', sub(pw(p2X, 2), pw(p2Y, 2)), 'sym4'); F_('QX^2+QY^2', add(pw(QX, 2), pw(QY, 2)), 'sym4')
F_('DX*e1X', mul(DX, e1X), 'disc4'); F_('DY*e1Y', mul(DY, e1Y), 'disc4'); F_('DX*e1Y', mul(DX, e1Y), 'disc4'); F_('DY*e1X', mul(DY, e1X), 'disc4')
F_('e3X*e1X', mul(e3X, e1X), 'sym4'); F_('e3Y*e1Y', mul(e3Y, e1Y), 'sym4'); F_('e3X*e1Y', mul(e3X, e1Y), 'sym4'); F_('e3Y*e1X', mul(e3Y, e1X), 'sym4')
F_('p3Y*e1Y', mul(p3Y, e1Y), 'sym4'); F_('p3X*e1X', mul(p3X, e1X), 'sym4'); F_('p3X*e1Y', mul(p3X, e1Y), 'sym4'); F_('p3Y*e1X', mul(p3Y, e1X), 'sym4')
F_('X0^4+X1^4+Y0^4+Y1^4', add(add(pw(X[0],4), pw(X[1],4)), add(pw(Y[0],4), pw(Y[1],4))), 'K3cone')
F_('X0^4+X1^4+X2^4+Y0^4', add(add(pw(X[0],4), pw(X[1],4)), add(pw(X[2],4), pw(Y[0],4))), 'K3cone')
F_('X0^4+Y0^4+Y1^4+Y2^4', add(pw(X[0],4), p4Y), 'K3cone')
F_('(X0^2-Y0^2)(X1^2-Y1^2)', mul(sub(pw(X[0],2), pw(Y[0],2)), sub(pw(X[1],2), pw(Y[1],2))), 'quad4')
F_('QX^2+QX*QY+QY^2', add(add(pw(QX,2), mul(QX,QY)), pw(QY,2)), 'sym4')
F_('(QX-QY)(QX+2QY)', mul(sub(QX,QY), add(QX, sc(QY,2))), 'sym4')
F_('QX(QX+2QY)', mul(QX, add(QX, sc(QY,2))), 'sym4'); F_('QY(QX+2QY)', mul(QY, add(QX, sc(QY,2))), 'sym4')
F_('e1X^4+e1Y^4', add(pw(e1X,4), pw(e1Y,4)), 'sym4')
if __name__ == '__main__':
    sel = sys.argv[1] if len(sys.argv) > 1 else 'all'
    for label, D, fam in fams:
        d = core.degree(D)
        if sel == 'low' and d > 4: continue
        t0 = time.time()
        r = test_D(F, A, q, D, fam, note=label, fast=(d >= 5), label=label)
        print(f"{label:28s} d={d} e1={r['e1']} cons={r.get('consistent')} rank={r.get('rank')} null={r.get('nullity')} "
              f"rows={r.get('rows')} cols={r.get('cols')} ({time.time()-t0:.1f}s)", flush=True)
        if r.get('consistent'):
            print('!!!!! CONSISTENT', label, flush=True)
            if 'N' in r:
                json.dump({'D': core.to_str(D), 'N': [core.to_str(n) for n in r['N']], 'Q': core.to_str(r['Q'])},
                          open(f'{RES}/witness-candidate-{label}.json', 'w'), indent=1)
