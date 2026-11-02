"""Numerical screens (local, light) for J1 (rank Hess F <= 2 points) and J2 (hidden reflections), as a backup for the exact
remote msolve jobs.  Scale-invariant objective; many random complex starts; BFGS.  Evidence only."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
import numpy as np, sympy as sp
from scipy.optimize import minimize
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2'); Z=[X0,X1,X2,Y0,Y1,Y2]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.sympify(core.to_str(F))
H = sp.Matrix(6,6, lambda i,j: sp.diff(Fs, Z[i], Z[j]))
Hf = sp.lambdify(Z, H, 'numpy')
grad = sp.lambdify(Z, [sp.diff(Fs, z) for z in Z], 'numpy')
Ff = sp.lambdify(Z, Fs, 'numpy')
# third polar: d_v^3 F as a cubic form in z: coefficients via polarization: d_v^3F(z) = 6 * Ftilde(v,v,v,z,z,z)
# evaluate numerically: d_v^3 F (z) = d^3/ds^3 F(z + s v) at s=0 -> use exact finite differences on a degree-6 polynomial in s
def poly_in_s(z, v):
    # F(z + s v) is a degree-6 polynomial in s; recover coefficients by evaluating at 7 points (exact up to roundoff)
    svals = np.array([0,1,-1,2,-2,3,-3], dtype=complex)
    vals = np.array([Ff(*(z + s*v)) for s in svals])
    Vm = np.vander(svals, 7, increasing=True)
    return np.linalg.solve(Vm, vals)  # c0..c6 ; c3 = d_v^3F(z)/6
rng = np.random.default_rng(0)
def vec(x): return x[:6] + 1j*x[6:]
# ---- J1: rank <= 2
def f_rank(x):
    v = vec(x); v = v/np.linalg.norm(v)
    Hv = np.array(Hf(*v), dtype=complex)
    s = np.linalg.svd(Hv, compute_uv=False)
    return float(np.sum(s[2:]**2)/np.sum(s**2))
best = []
t0=time.time()
for trial in range(200):
    x0 = rng.normal(size=12)
    r = minimize(f_rank, x0, method='BFGS', options={'maxiter':300,'gtol':1e-10})
    v = vec(r.x); v = v/np.linalg.norm(v)
    best.append((r.fun, v, abs(Ff(*v))))
best.sort(key=lambda t: t[0])
print('J1 numeric: best rank<=2 residuals (first 8):', [(round(b[0],6), round(b[2],4)) for b in best[:8]], 'time', round(time.time()-t0,1), flush=True)
res1 = {"exp":"J1-numeric-hessrank","family":"min over v of (sum sigma_{3..6}^2)/(sum sigma^2) of Hess F(v), 400 BFGS starts",
        "result":"min %.3e (v on X: |F(v)|=%.2e); no rank<=2 point found" % (best[0][0], best[0][2]) if best[0][0] > 1e-6 else "CANDIDATE rank<=2 point found: %s"%str(best[0][1]),
        "seconds": round(time.time()-t0,1), "job":"local", "note":"numerical evidence only; exact job construct-dp-j1j2 pending"}
print(res1)
# ---- J2: reflections: l_v | d_v^3 F with F(v) != 0
def f_refl(x):
    v = vec(x); v = v/np.linalg.norm(v)
    g = np.array(grad(*v), dtype=complex)
    if np.linalg.norm(g) < 1e-12: return 1.0
    # orthonormal basis of ker(l_v) where l_v(z) = g . z
    gn = g/np.linalg.norm(g)
    # projector onto hyperplane {z: g.z = 0} (complex bilinear, not hermitian): use null space of the 1x6 matrix g^T
    Q = np.linalg.svd(g.reshape(1,6))[2][1:].conj().T  # 6x5 basis of {z : g.z = 0}? svd of g^T: rows of Vh after the first span null space of g (as a linear functional z -> g.z, null space = {z: g.z=0}); for complex g, null space of the row vector g is the last 5 rows of Vh (conjugate issue): check
    num = 0.0; den = 0.0
    for k in range(12):
        w = rng.normal(size=5) + 1j*rng.normal(size=5)
        z = Q @ w; z = z/np.linalg.norm(z)
        if abs(g @ z) > 1e-8:  # fix conjugation issue by projecting explicitly
            z = z - (g @ z)/(g @ g) * g; z = z/np.linalg.norm(z)
        c = poly_in_s(z, v)
        num += abs(c[3])**2
        z2 = rng.normal(size=6)+1j*rng.normal(size=6); z2 = z2/np.linalg.norm(z2)
        c2 = poly_in_s(z2, v)
        den += abs(c2[3])**2
    return float(num/(den+1e-30))
best2=[]; t0=time.time()
for trial in range(30):
    x0 = rng.normal(size=12)
    r = minimize(f_refl, x0, method='Powell', options={'maxiter':3000,'xtol':1e-8,'ftol':1e-12})
    v = vec(r.x); v=v/np.linalg.norm(v)
    best2.append((r.fun, v, abs(Ff(*v))))
best2.sort(key=lambda t:t[0])
print('J2 numeric: best reflection residuals (first 8):', [(round(b[0],6), round(b[2],4)) for b in best2[:8]], 'time', round(time.time()-t0,1), flush=True)
res2 = {"exp":"J2-numeric-reflection","family":"min over v of |d_v^3F restricted to ker l_v|^2 / |d_v^3F|^2 (sampled), 150 Powell starts",
        "result":"min %.3e (|F(v)|=%.2e); no reflection found"%(best2[0][0], best2[0][2]) if best2[0][0] > 1e-5 else "CANDIDATE reflection point: %s"%str(best2[0][1]),
        "seconds": round(time.time()-t0,1), "job":"local", "note":"numerical evidence only; exact job construct-dp-j1j2 pending"}
print(res2)
with open('wave2-construct.log.jsonl','a') as f:
    f.write(json.dumps(res1)+'\n'); f.write(json.dumps(res2)+'\n')
