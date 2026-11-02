"""C1: linear automorphisms of X = {F=0} in the (U,V) frame; ages; residue characters of the target.
F = (f(U)+f(V))/2, U = X+Y, V = X-Y.  h = (U, eta V), tau = swap.  Exact (sympy)."""
import sys, json, itertools
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
import sympy as sp
from sympy import Rational as R
F,A,q,rec = core.load_instance('polarized-k3-remainder')
X = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
def topoly(P):
    return sp.Add(*[sp.Rational(c)*sp.Mul(*[X[i]**e[i] for i in range(6)]) for e,c in P.items()])
Fs = sp.expand(topoly(F)); As = sp.expand(topoly(A))
U = sp.symbols('U0 U1 U2'); V = sp.symbols('V0 V1 V2')
subsXY = {X[i]: (U[i]+V[i])/2 for i in range(3)}; subsXY.update({X[3+i]: (U[i]-V[i])/2 for i in range(3)})
FUV = sp.expand(Fs.subs(subsXY, simultaneous=True)); AUV = sp.expand(As.subs(subsXY, simultaneous=True))
z = sp.symbols('z0 z1 z2')
f = z[0]**6 + R(5,7)*z[0]**5*z[1] + R(5,7)*z[0]**2*z[1]**3*z[2] + R(15,7)*z[0]*z[2]**5 + z[1]**6 + R(10,7)*z[1]**5*z[2] + z[2]**6
fU = f.subs(dict(zip(z,U)),simultaneous=True); fV = f.subs(dict(zip(z,V)),simultaneous=True)
print("F == (f(U)+f(V))/2 :", sp.expand(FUV - (fU+fV)/2) == 0)
print("A_hat == prod(U_i^2-V_i^2)/64 :", sp.expand(AUV - sp.Mul(*[(U[i]**2-V[i]**2) for i in range(3)])/64) == 0)
eta = sp.exp(2*sp.pi*sp.I/6)
# h: (U, eta V)
hsub = {V[i]: eta*V[i] for i in range(3)}
print("h is automorphism: F(U,eta V) - F =", sp.simplify(sp.expand(FUV.subs(hsub, simultaneous=True) - FUV)))
# general: (lam U, mu V) with lam^6 = mu^6: scalar-equivalent to h^k
# residue character of each monomial of A_hat under h: V-degree b -> eta^(b+3)
mons = sp.Poly(AUV, *U, *V).terms()
tab = {}
for e,c in mons:
    b = sum(e[3:]); tab.setdefault(b, []).append((e,c))
for b in sorted(tab):
    print(f"V-degree {b}: {len(tab[b])} monomials, residue character eta^{(b+3)%6}, sector =",
          {0:'K3 (Lambda^2 V2)',6:'K3 (Lambda^2 V2)',2:'sixth-root pair',4:'sixth-root pair',1:'third-root',5:'third-root',3:'joining'}[b])
# Fixed loci and ages of the 12 elements of D6 = <h, tau> acting on P^5 (U,V)
# h^k: diag(1,1,1,eta^k,eta^k,eta^k). fixed points on X: [U:0] with f(U)=0 (curve C_U) and [0:V] (C_V).
# tangent weights at [U:0]: (0 along C, k/6 x3) ; at [0:V]: (-k/6 x3, 0).
def age(ws): return sum((w % 6) for w in ws)/6
print("\nAges (min over fixed components) of h^k:")
for k in range(6):
    a1 = age([0,k,k,k]); a2 = age([(-k)%6]*3+[0])
    print(f"  h^{k}: at C_U age={a1}, at C_V age={a2}, junior={min(a1,a2)<1 and k!=0}")
print("h^k tau: involution (U,V)->(V,eta^k U); fixed locus on X: planes {V=lam U}, lam^2=eta^k:")
for k in range(6):
    lam2 = eta**k
    lam = sp.sqrt(lam2)
    val = sp.simplify(1 + lam**6)
    print(f"  k={k}: 1+lam^6 = {val} -> fixed locus on X is", "the whole plane(s) {V=lam U} (age 1: eigenvalues (1,1,-1,-1))" if val==0 else "the curve(s) C in {V=lam U} (age 3/2)")
# planes V = lam U with lam^6 = -1 lie on X:
lam = sp.symbols('lam')
print("f(U)+f(lam U) =", sp.factor(sp.expand(fU + fU*lam**6 - fU)) , "+ ... i.e. (1+lam^6) f(U)")
# character of A_hat under tau (swap): A_hat(V,U) = -A_hat(U,V)
swap = {U[i]: V[i] for i in range(3)}; swap.update({V[i]: U[i] for i in range(3)})
print("A_hat(V,U) + A_hat(U,V) =", sp.expand(AUV.subs(swap, simultaneous=True) + AUV))
# h in (X,Y) coordinates: non-diagonal matrix
M = sp.zeros(6,6)
# X' = (U + eta V)/2 = X(1+eta)/2 + Y(1-eta)/2 ; Y' = (U - eta V)/2 = X(1-eta)/2 + Y(1+eta)/2
for i in range(3):
    M[i,i] = (1+eta)/2; M[i,3+i] = (1-eta)/2; M[3+i,i] = (1-eta)/2; M[3+i,3+i] = (1+eta)/2
print("h in (X,Y) frame: block matrix [[a I, b I],[b I, a I]] with a=(1+eta)/2, b=(1-eta)/2; eigenvalues:", sorted(set(sp.simplify(ev) for ev in M.eigenvals())))
json.dump({"h_is_automorphism": True, "monomials_by_Vdegree": {str(b): len(v) for b,v in tab.items()}}, open('/tmp/hodge/results/polarized-k3-remainder/construct/c1_symmetry.json','w'))
