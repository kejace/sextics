"""Enumerate the 405 planes of the Fermat cubic C: sum U_i^3 = 0 (matching of {0..5} into 3
pairs x three cube roots of unity per pair), their pullbacks S_Pi = phi^{-1}(Pi) (K3 (2,2,2)
complete intersections in P^5, all tau-stable), the plane incidence (line/point/empty) and whether
tau acts freely on S_Pi (no common zero of the three restricted conics on Y=0 and on X=0).
Works mod p = 2^31-1 (p = 1 mod 3, zeta = cube root of unity).  Output: w2c_planes.json
"""
import json, itertools, sys
from flint import nmod_mat
P = 2147483647
PRIMES = [7, 11, 13, 17, 19, 23]
CROSS = [((1,2),(3,4)), ((0,2),(4,5)), ((0,1),(3,5)), ((0,1),(4,5)), ((1,2),(3,5)), ((0,2),(3,4))]
def inv(a): return pow(a, P-2, P)
# cube root of unity mod p
ZETA = None
for g in range(2, 1000):
    z = pow(g, (P-1)//3, P)
    if z != 1: ZETA = z; break
assert (ZETA*ZETA*ZETA) % P == 1 and ZETA != 1
def quadrics():
    Q = []
    for i in range(6):
        c = (5 * inv(PRIMES[i])) % P
        q = {}
        e = [0]*6; e[i] = 2; q[tuple(e)] = 1
        for (a,b) in CROSS[i]:
            e = [0]*6; e[a] += 1; e[b] += 1; q[tuple(e)] = (q.get(tuple(e),0) + c) % P
        Q.append(q)
    return Q
QS = quadrics()
def matchings(s):
    if not s: yield []; return
    a = s[0]
    for i in range(1, len(s)):
        b = s[i]; rest = s[1:i] + s[i+1:]
        for m in matchings(rest): yield [(a,b)] + m
def planes():
    out = []
    for m in matchings(list(range(6))):
        for roots in itertools.product(range(3), repeat=3):
            # linear forms U_a + zeta^r U_b
            L = []
            for (a,b), r in zip(m, roots):
                row = [0]*6; row[a] = 1; row[b] = pow(ZETA, r, P); L.append(row)
            out.append({"matching": m, "roots": list(roots), "L": L})
    return out
PL = planes()
assert len(PL) == 405
def plane_dim_intersection(L1, L2):
    M = nmod_mat(L1 + L2, P)
    return 6 - M.rank() - 1   # projective dimension of the intersection (-1 = empty)
def restricted_conics_regular(L, block):
    # block: indices kept (X: 0,1,2 ; Y: 3,4,5); other coords set to 0; test whether the 3 conics
    # q_k = sum_j L[k][j] Q_j restricted to the block have no common zero: rank of q_k*S_2 on S_4 = 15
    keep = block
    def restrict(q):
        r = {}
        for e, v in q.items():
            if all(e[j] == 0 for j in range(6) if j not in keep):
                r[tuple(e[j] for j in keep)] = v
        return r
    qs = []
    for row in L:
        q = {}
        for j, c in enumerate(row):
            if c:
                for e, v in QS[j].items(): q[e] = (q.get(e, 0) + c*v) % P
        qs.append(restrict(q))
    def monos3(d):
        return [(a, b, d-a-b) for a in range(d+1) for b in range(d+1-a)]
    m2 = monos3(2); m4 = monos3(4); idx = {m:i for i,m in enumerate(m4)}
    rows = []
    for q in qs:
        for m in m2:
            r = [0]*15
            for e, v in q.items():
                k = tuple(a+b for a,b in zip(e, m)); r[idx[k]] = (r[idx[k]] + v) % P
            rows.append(r)
    return nmod_mat(rows, P).rank() == 15
def main():
    res = {"zeta": ZETA, "p": P, "planes": PL}
    inc = {}
    counts = {}
    for i in range(405):
        for j in range(i+1, 405):
            d = plane_dim_intersection(PL[i]["L"], PL[j]["L"])
            counts[d] = counts.get(d, 0) + 1
            if d >= 0: inc[f"{i},{j}"] = d
    res["incidence_counts"] = counts
    res["incidence"] = inc
    free = []
    for i, pl in enumerate(PL):
        fx = restricted_conics_regular(pl["L"], [0,1,2])
        fy = restricted_conics_regular(pl["L"], [3,4,5])
        free.append([fx, fy])
    res["tau_free_on_K3"] = free
    res["n_tau_free"] = sum(1 for a,b in free if a and b)
    json.dump(res, open("/tmp/hodge/results/quadratic-cubic-pullback/w2c_planes.json", "w"))
    print("incidence counts (proj dim of Pi cap Pi'):", counts)
    print("tau-free K3s:", res["n_tau_free"], "of 405")
if __name__ == "__main__": main()
