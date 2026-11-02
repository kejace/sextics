"""Full Fano variety of lines of W in a random chart of G(2,5): random T in GL5(F_p) applied to (x0,x1,y0,y1,y2),
line (s,t) -> T^{-1}(s, t, a2 s+b2 t, a3 s+b3 t, a4 s+b4 t); unknowns a2,a3,a4,b2,b3,b4,c0..c3 (10), 7 equations."""
import sys, random, sympy as sp
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave4')
from lines import *
p = int(sys.argv[1]); seed = int(sys.argv[2]); random.seed(seed)
a2,a3,a4,b2,b3,b4 = sp.symbols('a2 a3 a4 b2 b3 b4'); cc = sp.symbols('c0 c1 c2 c3')
cub = cc[0]*s**3 + cc[1]*s**2*t + cc[2]*s*t**2 + cc[3]*t**3
while True:
    T = sp.Matrix(5, 5, lambda i, j: random.randrange(-5, 6))
    if T.det() != 0: break
w = sp.Matrix([s, t, a2*s+b2*t, a3*s+b3*t, a4*s+b4*t])
lin = list(T * w)          # new coords = T * (chart coords)  (any invertible T gives a random chart)
eqs = incidence(lin, cub)
vars_ = [a2,a3,a4,b2,b3,b4] + list(cc)
fn = f'fano_full_{p}_{seed}'
msolve_file(fn + '.ms', eqs, vars_, p)
singular_file(fn + '.sing', eqs, vars_, p,
  'option(redSB);\nideal J = std(I);\nprint("DIM"); print(dim(J)); print("MULT"); print(mult(J));\n'
  'write(":w ' + fn + '_gb.txt", J);\n'
  'LIB "primdec.lib";\nlist L = minAssGTZ(J);\nprint("NCOMP"); print(size(L));\n'
  'int i; for (i=1; i<=size(L); i++) { ideal K = std(L[i]); print("COMP"); print(i); print(dim(K)); print(mult(K)); print(size(L[i])); write(":w ' + fn + '_comp"+string(i)+".txt", L[i]); }\nquit;\n')
print('wrote', fn, 'T=', T.tolist())
open(fn + '_T.txt','w').write(str(T.tolist()))
