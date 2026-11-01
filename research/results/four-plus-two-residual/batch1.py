"""Batch 1: structural denominators (d<=6)."""
from h import *
Fs = core.to_str
X0,X1,X2,Y0,Y1,Y2 = [P(v) for v in core.VARS]
H = P('X0**6 + 5*X0**5*X1/7 + 5*X0**3*X1*Y0*Y1/7 + 10*X0*X1**5/7 + 10*X0*X1*Y0**3*Y1/7 + X1**6 + Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6')
B = P('X2**6 + 15*X2**4*Y2**2/7 + 10*X2**2*Y2**4/7 + Y2**6')
f = P('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6')
g = P('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6')
K = P('5*X0**3*X1*Y0*Y1/7 + 10*X0*X1*Y0**3*Y1/7')
assert core.sub(core.add(core.add(f,g),core.add(K,B)), F) == {}
tests = [
 ('X0**2 + 2*Y0**2', 'structural', 'quadratic factor of coupling'),
 ('X0**2 - 2*Y0**2', 'structural', 'conjugate quadratic'),
 ('X0**2 + Y0**2', 'structural', ''),
 ('X0*X1 + Y0*Y1', 'structural', ''),
 ('X0*X1*Y0*Y1', 'structural', 'coupling monomial part'),
 (Fs(K), 'structural', 'coupling term K'),
 (Fs(B), 'structural', 'binary block B'),
 (Fs(H), 'structural', 'core block H'),
 (Fs(f), 'structural', 'f(X0,X1)'),
 (Fs(g), 'structural', 'g(Y0,Y1)'),
 (Fs(core.add(f,g)), 'structural', 'f+g (H minus coupling)'),
 (Fs(core.add(core.add(f,g),B)), 'structural', 'F minus coupling = three binary blocks'),
 ('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6', 'structural', 'pure power part'),
 (Fs(core.add(H, core.scale(B, -1))), 'structural', 'H - B'),
 (Fs(core.add(f, B)), 'structural', 'f + B'),
 (Fs(core.add(g, B)), 'structural', 'g + B'),
 (Fs(core.add(f, g)), 'structural', 'f+g again? skip', ),
]
# partial derivatives (pass E1 trivially)
for i, v in enumerate(core.VARS):
    tests.append((Fs(core.pderiv(F, i)), 'partials', f'dF/d{v}'))
# block radial derivatives
for (i,j,name) in [(0,3,'X0,Y0'),(1,4,'X1,Y1'),(2,5,'X2,Y2'),(0,1,'X0,X1'),(3,4,'Y0,Y1')]:
    R = core.add(core.mul(P(core.VARS[i]), core.pderiv(F,i)), core.mul(P(core.VARS[j]), core.pderiv(F,j)))
    tests.append((Fs(R), 'radial', f'({name}) radial derivative of F'))
seen = set()
for expr, fam, note in tests:
    if expr in seen: continue
    seen.add(expr)
    r = run(expr, fam, note=note)
    print({k: r[k] for k in ('D','d','e1','consistent','rank','nullity','seconds','note')}, flush=True)
