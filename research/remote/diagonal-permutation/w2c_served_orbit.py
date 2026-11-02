"""Helper for the search lane (construct lens, §0.9 of the report): given a denominator D, decide whether the S3-orbit product
prod_g D∘g serves A, i.e. whether A ∈ Q[S3]·served(D) + J_6^-.  Computes served(D) with fixed_d.target_image (mod p), applies the
six simultaneous permutations to its basis, and tests membership of A in the span modulo J(F)_6^- (Y-odd).
Usage: python3 w2c_served_orbit.py "<D expression>" [p]"""
import sys, json, time, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian
from hodge.jacobian import jacobian_generators, _vectorize, reduce_mod
from flint import nmod_mat
F, A, q, rec = core.load_instance('diagonal-permutation')
def permute(poly, perm):
    # simultaneous permutation of pairs: (X_i,Y_i) -> (X_perm[i], Y_perm[i])
    out = {}
    for e, c in poly.items():
        ne = [0]*6
        for i in range(3):
            ne[perm[i]] += e[i]; ne[3+perm[i]] += e[3+i]
        out[tuple(ne)] = (out.get(tuple(ne), 0) + c)
    return {e:c for e,c in out.items() if c}
def orbit_test(D, p=2147483647, verbose=True):
    t0 = time.time()
    rank, basis, S = fixed_d.target_image(F, q, D, p=p)
    gens = [g[2] for g in jacobian_generators(F, 6, p, 1)]
    vecs = list(gens)
    for b in basis:
        for perm in itertools.permutations(range(3)):
            vecs.append(reduce_mod(permute(b, perm), p))
    rows = _vectorize(vecs + [reduce_mod(A, p)])
    def mat(cols):
        R, C = len(rows), len(cols); flat=[0]*(R*C)
        for j,col in enumerate(cols):
            for e,v in col.items(): flat[rows[e]*C+j] = v % p
        return nmod_mat(R, C, flat, p)
    r1 = mat(vecs).rank(); r2 = mat(vecs + [reduce_mod(A, p)]).rank()
    rJ = mat(gens).rank()
    res = {"served_dim": rank, "served_primitive_dim": rank - rJ if rank >= rJ else rank, "orbit_span_primitive_dim": r1 - rJ,
           "A_in_orbit_span": (r2 == r1), "seconds": round(time.time()-t0,1), "prime": p}
    if verbose: print(res)
    return res
if __name__ == '__main__':
    Dexpr = sys.argv[1]; p = int(sys.argv[2]) if len(sys.argv) > 2 else 2147483647
    D = core.parse(Dexpr)
    res = orbit_test(D, p)
    res.update({"exp":"served-orbit-test","D":Dexpr,"d":core.degree(D),"job":"local"})
    with open('wave2-construct.log.jsonl','a') as f: f.write(json.dumps(res)+'\n')
