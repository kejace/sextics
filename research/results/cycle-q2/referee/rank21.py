"""Exact rank over Q of the 21 second partials F_ij of F (cycle-q2) -- referee check (no engine)."""
from fractions import Fraction
import itertools, ast
VARS = ['X0','X1','X2','Y0','Y1','Y2']
def parse(s):
    def ev(n):
        if isinstance(n, ast.Expression): return ev(n.body)
        if isinstance(n, ast.Constant): return {(0,)*6: Fraction(n.value)}
        if isinstance(n, ast.Name): e=[0]*6; e[VARS.index(n.id)]=1; return {tuple(e): Fraction(1)}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return {k:-v for k,v in ev(n.operand).items()}
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): r=dict(a); [r.__setitem__(k, r.get(k,0)+v) for k,v in b.items()]; return {k:v for k,v in r.items() if v}
            if isinstance(n.op, ast.Sub): r=dict(a); [r.__setitem__(k, r.get(k,0)-v) for k,v in b.items()]; return {k:v for k,v in r.items() if v}
            if isinstance(n.op, ast.Mult):
                r={}
                for k1,v1 in a.items():
                    for k2,v2 in b.items():
                        k=tuple(x+y for x,y in zip(k1,k2)); r[k]=r.get(k,0)+v1*v2
                return {k:v for k,v in r.items() if v}
            if isinstance(n.op, ast.Div):
                assert len(b)==1 and list(b)[0]==(0,)*6; c=list(b.values())[0]; return {k:v/c for k,v in a.items()}
            if isinstance(n.op, ast.Pow):
                e=int(list(b.values())[0]); r={(0,)*6:Fraction(1)}
                for _ in range(e): r=ev(ast.BinOp(left=ast.Constant(0),op=ast.Mult(),right=ast.Constant(0))) if False else mulp(r,a)
                return r
        raise ValueError(ast.dump(n))
    def mulp(a,b):
        r={}
        for k1,v1 in a.items():
            for k2,v2 in b.items():
                k=tuple(x+y for x,y in zip(k1,k2)); r[k]=r.get(k,0)+v1*v2
        return {k:v for k,v in r.items() if v}
    return ev(ast.parse(s, mode='eval'))
def pderiv(P,i):
    r={}
    for k,v in P.items():
        if k[i]:
            e=list(k); e[i]-=1; r[tuple(e)]=r.get(tuple(e),0)+v*k[i]
    return r
F = parse(open('/home/kejace/devel/sextics/data/cycle-q2/instance.rec').read().split('F: ')[1].split('\n')[0])
Fij = [pderiv(pderiv(F,i),j) for i in range(6) for j in range(i,6)]
mons = sorted({e for P in Fij for e in P}); idx={e:t for t,e in enumerate(mons)}
# exact rank by fraction-free Gaussian elimination
M = [[Fraction(0)]*len(mons) for _ in Fij]
for r,P in enumerate(Fij):
    for e,c in P.items(): M[r][idx[e]] = c
rank=0; rows=len(M); cols=len(mons)
for c in range(cols):
    piv = next((r for r in range(rank, rows) if M[r][c] != 0), None)
    if piv is None: continue
    M[rank], M[piv] = M[piv], M[rank]
    for r in range(rows):
        if r != rank and M[r][c] != 0:
            f = M[r][c]/M[rank][c]; M[r] = [a - f*b for a,b in zip(M[r], M[rank])]
    rank += 1
print('F terms', len(F), '| 21 second partials span rank', rank, 'in', cols, 'quartic monomials')
