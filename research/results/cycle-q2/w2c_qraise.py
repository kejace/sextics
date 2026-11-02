#!/usr/bin/env python3
"""(iv) Verify the q-raising identity symbolically (dossier §5 P2 / A.7), in the contract's cleared form.
Setting: F_s = F + s*M (M any polynomial), B with A = M*B.  Suppose (D_s, n_s, Q_s) is a q=1 certificate
for (F_s, B) with coefficients polynomial in s:  B D^2 - D n(F_s) + F_s (D div n - n(D)) = F_s^2 Q   (*)
where n(F_s) = sum n_i d_i F_s.  Claim: with v = n/D the vector field, w := (M v - F_s v_dot)/2 (dot = d/ds)
is a q=2 certificate for (F_s, A): A - 2 w(F_s) + F_s div w = -(1/2) F_s^3 Q_dot (rational form), hence
with denominator D^2: numerators W_i = (M n_i D - F_s (n_i_dot D - n_i D_dot))/2 and quotient -Q_dot/2.
We check the cleared identity  A W^2... no: cleared form with denominator E = D^2:
   A E^2 - 2 E sum W_i F_i + F_s (E div W - W(E)) = F_s^3 * (-(Q_dot/2) * ???)
To avoid guessing the exact cleared normalisation we verify the rational identity directly with sympy on a
2-variable toy where (*) holds identically: we construct (*) by choosing n, D, F_s, B freely and DEFINING
Q := (B D^2 - D n(F_s) + F_s (D div n - n(D)))/F_s^2 -- this requires F_s^2 | numerator, which generic
choices violate; so instead we verify the derivative identity as a formal consequence:
   let E1(s) := B - v(F_s) + F_s div v  (rational in z, s), with v = n/D.  Then
   d/ds E1 = -v(M) - v_dot(F_s) + M div v + F_s div v_dot,
   and  A - 2 w(F_s) + F_s div w  with w = (M v - F_s v_dot)/2  equals  M*E1 - (F_s/2) d/ds E1  identically.
Hence if E1 = F_s^2 Q then A - 2w(F_s) + F_s div w = F_s^2 (M Q - F_s Q_dot - M Q ... ) -- we let sympy do it."""
import sympy as sp
x,y,s=sp.symbols('x y s')
Z=[x,y]
def rnd(deg,seed):
    import random; random.seed(seed)
    return sum(random.randint(-3,3)*x**i*y**(deg-i) for i in range(deg+1))
F0=rnd(4,1); M=rnd(4,2); Fs=F0+s*M
B=rnd(2,3); A=M*B
D=rnd(2,4)+s*rnd(2,5)
n=[rnd(3,6)+s*rnd(3,7), rnd(3,8)+s*rnd(3,9)]
v=[ni/D for ni in n]
def vec_apply(vv,P): return sum(vi*sp.diff(P,zi) for vi,zi in zip(vv,Z))
def div(vv): return sum(sp.diff(vi,zi) for vi,zi in zip(vv,Z))
E1 = B - vec_apply(v,Fs) + Fs*div(v)
vdot=[sp.diff(vi,s) for vi in v]
w=[(M*vi - Fs*vdi)/2 for vi,vdi in zip(v,vdot)]
E2 = A - 2*vec_apply(w,Fs) + Fs*div(w)
claim = sp.simplify(E2 - (M*E1 - Fs*sp.diff(E1,s)/2))
print("E2 - (M*E1 - (F_s/2) dE1/ds) =", claim)
# consequence: if E1 = F_s^2 Q then E2 = M F_s^2 Q - (F_s/2)(2 F_s M Q + F_s^2 Qdot) = -(1/2) F_s^3 Qdot
Q=sp.Function('Q')(x,y,s)
E1q=Fs**2*Q
print("if E1=F_s^2 Q:", sp.simplify(M*E1q - Fs*sp.diff(E1q,s)/2 + sp.Rational(1,2)*Fs**3*sp.diff(Q,s)))
