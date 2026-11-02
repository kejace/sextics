# Generate the 13x13 profile-A incidence system on W_{23,f} (H = 7f(x,a,y) + 7S^2 + 7T^3, x = 4t) as msolve input.
# usage: gen_incidence_msolve.py <k0> <k1> <k2> <k3> <k4> <char> <outfile>
import sys, sympy as sp
k0,k1,k2,k3,k4,ch = [int(v) for v in sys.argv[1:7]]; out=sys.argv[7]
t=sp.symbols('t'); c=sp.symbols('c0:13')
x=4*t; a=c[0]+c[1]*t+c[2]*t**2; y=k0+c[3]*t+c[4]*t**2
S=k1+c[5]*t+c[6]*t**2+c[7]*t**3+c[8]*t**4+c[9]*t**5+k3*t**6
T=k2+c[10]*t+c[11]*t**2+c[12]*t**3+k4*t**4
f7 = 7*x**6 + 5*x**5*a + 5*x**2*a**3*y + 15*x*y**5 + 7*a**6 + 10*a**5*y + 7*y**6   # 7 f(x,a,y)
H = sp.expand(f7 + 7*S**2 + 7*T**3)
P = sp.Poly(H, t)
eqs=[sp.expand(P.coeff_monomial(t**i)) for i in range(13)]
assert P.degree()==12
def fmt(e):
    s=sp.sstr(e).replace('**','^').replace(' ','')
    return s
with open(out,'w') as fo:
    fo.write(','.join(f'c{i}' for i in range(13))+'\n'); fo.write(f'{ch}\n')
    fo.write(',\n'.join(fmt(e) for e in eqs)+'\n')
print('written', out, 'degrees', [sp.Poly(e,*c).total_degree() for e in eqs])
