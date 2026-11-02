# Tangent space of the Fano scheme of lines F_1(X) at the Fermat lines lying on the ten Fermat surfaces X ∩ {Z_i=Z_j=0}.
# Exact arithmetic mod p = 65521 (p = 1 mod 12, so zeta with zeta^6 = -1 exists; 7 invertible).
import sys, json, itertools, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
p = 65521
F, A, q, rec = core.load_instance('reflection-coupled-residual')
# F as dict {exp: Fraction}; reduce mod p
def modp(fr): return (fr.numerator % p) * pow(fr.denominator % p, p-2, p) % p
Fm = {e: modp(c) for e, c in F.items()}
def pderiv(P, i):
    out = {}
    for e, c in P.items():
        if e[i] > 0:
            e2 = list(e); e2[i] -= 1
            out[tuple(e2)] = (out.get(tuple(e2), 0) + c * e[i]) % p
    return out
Fi = [pderiv(Fm, i) for i in range(6)]
def polymul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % p
    return r
def eval_line(P, pv, qv):
    # returns coefficients in t of P(pv + t qv)
    tot = [0]*7
    for e, c in P.items():
        term = [c]
        for j in range(6):
            lin = [pv[j] % p, qv[j] % p]
            for _ in range(e[j]): term = polymul(term, lin)
        for k, x in enumerate(term): tot[k] = (tot[k] + x) % p
    return tot
def rank_mod_p(M):
    M = [row[:] for row in M]; r = 0; rows = len(M); cols = len(M[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [x*inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(x - f*y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r
def tangent_dim(pv, qv):
    # chart: pick indices r,s with invertible minor; unknowns = other coordinates of p and q (after normalisation)
    for r, s in itertools.combinations(range(6), 2):
        det = (pv[r]*qv[s] - pv[s]*qv[r]) % p
        if det: break
    # normalise: new p' = (qs p - ps q)/det, q' = (-qr p + pr q)/det  so that p'_r=1,p'_s=0,q'_r=0,q'_s=1
    inv = pow(det, p-2, p)
    pn = [((qv[s]*pv[j] - pv[s]*qv[j]) * inv) % p for j in range(6)]
    qn = [((-qv[r]*pv[j] + pv[r]*qv[j]) * inv) % p for j in range(6)]
    assert pn[r] == 1 and pn[s] == 0 and qn[r] == 0 and qn[s] == 1
    assert all(x == 0 for x in eval_line(Fm, pn, qn)), 'line not on X'
    free = [j for j in range(6) if j not in (r, s)]
    # Jacobian: rows k=0..6, columns: dp_j (j in free), dq_j (j in free)
    J = [[0]*8 for _ in range(7)]
    for ci, j in enumerate(free):
        fj = eval_line(Fi[j], pn, qn)  # degree <= 5 -> entries 0..5
        for k in range(7):
            J[k][ci] = fj[k] if k < len(fj) else 0          # d/dp_j of [t^k]
            J[k][4+ci] = fj[k-1] if 1 <= k <= len(fj) else 0  # d/dq_j of [t^k] = [t^{k-1}] F_j
    return 8 - rank_mod_p(J)
# sixth roots of -1 mod p
roots = [z for z in range(1, p) if pow(z, 6, p) == p-1]
assert len(roots) == 6
results = {}
t0 = time.time()
lines_seen = set()
for blk_pair in [(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(0,4),(1,3),(1,4)]:
    rest = [j for j in range(6) if j not in blk_pair]
    dims = {}
    for pairing in [((rest[0],rest[1]),(rest[2],rest[3])), ((rest[0],rest[2]),(rest[1],rest[3])), ((rest[0],rest[3]),(rest[1],rest[2]))]:
        (k,l),(m,n) = pairing
        for z1 in roots:
            for z2 in roots:
                pv = [0]*6; qv = [0]*6
                pv[k] = z1; pv[l] = 1; qv[m] = z2; qv[n] = 1
                key = (tuple(pv), tuple(qv))
                d = tangent_dim(pv, qv)
                dims[d] = dims.get(d, 0) + 1
                lines_seen.add(key)
    results[str(blk_pair)] = dims
    print(blk_pair, 'tangent-space dims of F_1(X) at the 108 lines:', dims, flush=True)
print('total lines', len(lines_seen), 'time', round(time.time()-t0,1))
json.dump({'p': p, 'tangent_dims_by_surface': results, 'n_lines': len(lines_seen)}, open('/tmp/hodge/results/reflection-coupled-residual/construct/lines_tangent.json','w'), indent=1)
