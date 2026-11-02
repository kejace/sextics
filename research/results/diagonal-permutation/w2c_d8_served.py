"""Construct lens: served space of the unique degree-8 E1-admissible denominator D8 (and its (01)-symmetrisation), via
fixed_d.target_image, plus the S3-orbit-span test for A (w2c_served_orbit.orbit_test).  Written for the remote (dense nullspace
of a ~40k x 7.7k matrix mod p).  Outputs: d8_served.out (this stdout), log lines in wave2-construct.log.jsonl."""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hodge import core, fixed_d, jacobian
import sympy as sp
from w2c_served_orbit import orbit_test, permute
F, A, q, rec = core.load_instance('diagonal-permutation')
p = 2147483647
D8 = core.parse(json.load(open('D8_e1kernel.json'))['D'])
D8s = core.scale(core.add(D8, permute(D8, (1,0,2))), sp.Rational(1,2)) if hasattr(core,'scale') else None
if D8s is None:
    from fractions import Fraction
    tmp = core.add(D8, permute(D8, (1,0,2))); D8s = {e: Fraction(c)/2 for e,c in tmp.items()}
log = open('wave2-construct.log.jsonl','a')
for name, D in [('D8 (rref representative, 414 terms)', D8), ('D8 (01)-symmetrised', D8s)]:
    t0 = time.time()
    e1 = fixed_d.e1_solvable(F, A, q, D, p) if 'e1_solvable' in dir(fixed_d) else None
    print(name, 'deg', core.degree(D), 'terms', len(D), 'E1 for A:', e1, flush=True)
    try:
        res = orbit_test(D, p)
    except Exception as ex:
        res = {"error": repr(ex)}
    res.update({"exp":"D8-served-space","D":name,"d":8,"prime":p,"e1_A":str(e1),"seconds":round(time.time()-t0,1),"job":"construct-dp-d8served",
                "note":"served_primitive_dim = number of independent tau-invariant primitive (3,1)-classes supported on V(D8); A_in_orbit_span = A served by the S3-orbit product"})
    print(res, flush=True); log.write(json.dumps(res)+"\n"); log.flush()
