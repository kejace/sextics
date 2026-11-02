"""Wave-2 construct lens, experiment X1: kernel of multiplication by candidate degree-6 'cone' denominators
on the Y-odd part R_6^- of the Jacobian ring (E1 for the DIVISOR: served targets must lie in this kernel).
Candidates: cones over the tau-fixed curves (sum X_i^6, sum Y_i^6), transposition cones F|_{r=t=0}, Fermat control."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, jacobian
import sympy as sp
X0,X1,X2,Y0,Y1,Y2 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
V = [X0,X1,X2,Y0,Y1,Y2]
F, A, q, rec = core.load_instance('diagonal-permutation')
Fs = sp.sympify(rec['F'].replace('\n',' ')) if isinstance(rec, dict) and 'F' in rec else None
if Fs is None:
    Fs = sp.sympify(core.to_str(F))
p = 2147483647
def tostr(e):
    return str(sp.expand(e))
def cone_transposition(i, j):
    Xi, Xj, Yi, Yj = V[i], V[j], V[3+i], V[3+j]
    sub = {Xi: (Xi+Xj)/2, Xj: (Xi+Xj)/2, Yi: (Yi+Yj)/2, Yj: (Yi+Yj)/2}
    return sp.expand(Fs.subs(sub, simultaneous=True))
cands = {
 'sumX6': X0**6+X1**6+X2**6,
 'sumY6': Y0**6+Y1**6+Y2**6,
 'sumX6+sumY6': X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6,
 'cone01': cone_transposition(0,1),
 'cone02': cone_transposition(0,2),
 'pair0': X0**6+Y0**6,
 'F42': sp.expand(Fs - (X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6)),  # coupling part, for comparison
}
out = []
log = open('/tmp/hodge/results/diagonal-permutation/wave2-construct.log.jsonl','a')
for name, D in cands.items():
    Dd = core.parse(tostr(D))
    t0 = time.time()
    r = jacobian.mult_kernel(F, Dd, 6, p, parity=1)
    dt = time.time()-t0
    rec_ = {"exp":"X1-divisorE1-kernel","D":name,"d":6,"prime":p,"result":"ker dim %d (ring_dim %d)"%(r['dim'], r['ring_dim']),
            "rank": r['ring_dim']-r['dim'],"seconds":round(dt,1),"job":"local","note":"kernel of [D] on R_6^- = upper bound for primitive targets served by V(D)"}
    print(rec_, flush=True)
    log.write(json.dumps(rec_)+"\n"); log.flush()
    out.append((name, r['dim']))
# Fermat control
Ff = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
r = jacobian.mult_kernel(Ff, core.parse('X0**6+X1**6+X2**6'), 6, p, parity=1)
rec_ = {"exp":"X1-divisorE1-kernel","D":"sumX6 on FERMAT (control)","d":6,"prime":p,"result":"ker dim %d (ring_dim %d)"%(r['dim'], r['ring_dim']),"job":"local","note":"control: sumX6 in J(Fermat) so kernel = everything"}
print(rec_, flush=True); log.write(json.dumps(rec_)+"\n"); log.flush()
