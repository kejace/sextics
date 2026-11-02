# Numerical check: local monodromy of the zeta_6-eigenpart L of H^1 of z^6 = (U-e1)(U-e2)(U-e3)
# around the collision e2 -> e1.  Forms w_j = U^j dU / z^5 = U^j prod (U-e_k)^(-5/6) dU, j=0,1.
import mpmath as mp, json
mp.mp.dps = 25
e1 = mp.mpf(0); e3 = mp.mpf(1); rho = mp.mpf('0.3')
def e2(phi): return e1 + rho*mp.expj(phi)
def parg(zv):  # principal arg
    return mp.arg(zv)
def I12(phi, j):
    # U = e1 + (e2-e1) u ; branches: (U-e1)^(-5/6) = rho^(-5/6) e^{-5i phi/6} u^(-5/6)
    # (U-e2)^(-5/6) = rho^(-5/6) e^{-5i phi/6} e^{-5 i pi/6} (1-u)^(-5/6) ; (U-e3): principal (stays near -1)
    d = e2(phi) - e1
    pref = d * rho**(mp.mpf(-5)/3) * mp.expj(-mp.mpf(5)*phi/3) * mp.expj(-5*mp.pi/6)
    def f(u):
        U = e1 + d*u
        return U**j * u**(mp.mpf(-5)/6) * (1-u)**(mp.mpf(-5)/6) * mp.exp(-mp.mpf(5)/6*(mp.log(abs(U-e3)) + 1j*parg(U-e3)))
    return pref * mp.quad(f, [0, 1])
def I31(phi, j):
    # path: e3 -> Q=(e1+3rho) along real axis ; arc |U-e1|=3rho from angle 0 to phi+pi ; radial to e1
    R3 = 3*rho; total = mp.mpc(0)
    # leg (i): from e3 to Q, U = x real in [R3,1]; args: (U-e1):0 ; (U-e2): principal ; (U-e3): pi
    def fi(x):
        U = x
        a1 = mp.log(x)                      # log(U-e1)
        a2 = mp.log(abs(U-e2(phi))) + 1j*parg(U-e2(phi))
        a3 = mp.log(abs(1-x)) + 1j*mp.pi
        return U**j * mp.exp(-mp.mpf(5)/6*(a1+a2+a3))
    total += -mp.quad(fi, mp.linspace(R3, 1, 12))          # integral from e3 to Q = -(integral from Q to e3) = -int_{R3}^{1}
    # leg (ii): arc U = e1 + R3 e^{i psi}, psi from 0 to phi+pi
    def fii(psi):
        U = e1 + R3*mp.expj(psi)
        a1 = mp.log(R3) + 1j*psi
        a2 = mp.log(rho) + 1j*psi + mp.log(abs(3 - mp.expj(phi-psi))) + 1j*parg(3 - mp.expj(phi-psi))
        # U-e2 = rho(3 e^{i psi} - e^{i phi}) = rho e^{i psi} (3 - e^{i(phi-psi)})
        w = 1 - R3*mp.expj(psi)             # U - e3 = -(1 - 3 rho e^{i psi})
        a3 = mp.log(abs(w)) + 1j*(mp.pi + parg(w))
        dU = 1j*R3*mp.expj(psi)
        return U**j * mp.exp(-mp.mpf(5)/6*(a1+a2+a3)) * dU
    total += mp.quad(fii, mp.linspace(0, phi+mp.pi, 48))
    # leg (iii): radial U = e1 + s e^{i(phi+pi)}, s from R3 to 0
    def fiii(s):
        U = e1 + s*mp.expj(phi+mp.pi)
        a1 = mp.log(s) + 1j*(phi+mp.pi)
        a2 = mp.log(s+rho) + 1j*(phi+mp.pi)   # U-e2 = -(s+rho) e^{i phi}
        w = 1 + s*mp.expj(phi)                 # U - e3 = -(1 + s e^{i phi})
        a3 = mp.log(abs(w)) + 1j*(mp.pi + parg(w))
        dU = mp.expj(phi+mp.pi)
        return U**j * mp.exp(-mp.mpf(5)/6*(a1+a2+a3)) * dU
    total += -mp.quad(fiii, mp.linspace(0, R3, 12))        # s runs from R3 down to 0
    return total
def Pi(phi):  # noqa
    return mp.matrix([[I12(phi,0), I12(phi,1)], [I31(phi,0), I31(phi,1)]])
P0 = Pi(mp.mpf(0)); P1 = Pi(2*mp.pi)
T = P1 * P0**-1
ev = mp.eig(T)[0]
zeta3 = mp.expj(2*mp.pi/3)
print("Pi(0) =", P0); print("T =", T)
print("eigenvalues:", ev)
print("|ev - 1|, |ev - zeta3|, |ev - zeta3^2| for each:", [[mp.nstr(abs(e-1),6), mp.nstr(abs(e-zeta3),6), mp.nstr(abs(e-zeta3**2),6)] for e in ev])
print("T^3 - I norm:", mp.nstr(mp.norm(T**3 - mp.eye(2)), 6), "  T^2 - I norm:", mp.nstr(mp.norm(T**2 - mp.eye(2)), 6))
print("det T:", T[0,0]*T[1,1]-T[0,1]*T[1,0])
json.dump({"rho": str(rho), "eigenvalues": [str(e) for e in ev], "T3_minus_I": str(mp.norm(T**3 - mp.eye(2))), "T2_minus_I": str(mp.norm(T**2-mp.eye(2))), "detT": str(T[0,0]*T[1,1]-T[0,1]*T[1,0])}, open("monodromy_check_fine.json","w"), indent=1)
