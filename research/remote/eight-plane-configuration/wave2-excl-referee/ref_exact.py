"""Exact (characteristic-0) referee checks, own code:
 (1) A notin J(F)_6: rank of the 18 odd generators Z_j F_i vs with A appended (fmpq).
 (2) rank of the 21 second partials F_ij (quartics) over Q: 21 => Hess_F(w,w)=0 forces w=0 (any parity), closes degree 5 given E1 at d=4.
 (3) Koszul-reduced E2 for D = dF/dX0 and D = dF/dY0 over Q, NO parity restriction (c_ij in S_1: 90 unknowns, Q in S_4: 126):
     rank M vs rank [M|b]  (fmpz after clearing the denominators of F; column scalings do not affect rank/consistency)."""
import sys, time
from fractions import Fraction
from math import lcm
from flint import fmpz_mat, fmpq_mat, fmpq
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl-referee')
import ref_lib as L

F, A, q = L.load()
Fi = [L.pderiv(F, i) for i in range(L.NV)]
Fij = [[L.pderiv(Fi[i], j) for j in range(L.NV)] for i in range(L.NV)]

def rank_exact(polys):
    cols = {}
    for P in polys:
        for e in P: cols.setdefault(e, len(cols))
    return L.fmpq_from_rows(polys, cols, len(cols)).rank()

# (1)
t0 = time.time()
gens = L.J_generators(Fi, 6, 1)
r1 = rank_exact(gens); r2 = rank_exact(gens + [A])
L.log({'check': 'exact: A notin J_6', 'rank_J6odd': r1, 'rank_J6odd_plus_A': r2, 'A_in_J': r2 == r1, 'seconds': round(time.time() - t0, 1),
       'note': 'S1 of the excl report; its text says dim J_6 = 126, the correct number is 36 (18 per parity); the argument (mod-5 Fermat + saturation) is fine'})
# (2)
t0 = time.time()
quads = [Fij[i][j] for i in range(L.NV) for j in range(i, L.NV)]
rk21 = rank_exact(quads)
rkXX = rank_exact([Fij[i][j] for i in range(3) for j in range(i, 3)])
rkYY = rank_exact([Fij[i][j] for i in range(3, 6) for j in range(i, 6)])
L.log({'check': 'exact: rank of second partials', 'rank_all_21': rk21, 'rank_FXX_6': rkXX, 'rank_FYY_6': rkYY, 'seconds': round(time.time() - t0, 1),
       'conclusion': 'Hess_F(w,w) = sum w_i w_j F_ij = 0 (as a quartic) forces all w_i w_j = 0, i.e. w = 0, for constant w of ANY parity mix' if rk21 == 21 else 'FAIL'})
# (3) Koszul-reduced E2 over Q (no parity restriction)
den = 1
for c in F.values(): den = lcm(den, Fraction(c).denominator)
Fz = {e: int(Fraction(c) * den) for e, c in F.items()}         # integer model L*F
Fzi = [L.pderiv(Fz, i) for i in range(L.NV)]
Fzij = [[L.pderiv(Fzi[i], j) for j in range(L.NV)] for i in range(L.NV)]
Az = {e: int(c) for e, c in A.items()}
def koszul_exact(k):
    """D = d(LF)/dZ_k, w = e_k. Columns: c_ij monomials (deg 1), Q monomials (deg 4). Rows: degree-10 monomials."""
    t0 = time.time()
    D = Fzi[k]; Di = [L.pderiv(D, i) for i in range(L.NV)]
    cols = []
    for i in range(L.NV):
        for j in range(i + 1, L.NV):
            for m in L.monomials(1):
                dm = [L.pderiv({m: 1}, t) for t in range(L.NV)]
                # D * (F_j d_i m - F_i d_j m) - m * (F_j D_i - F_i D_j)
                t1 = L.add(L.mul(Fzi[j], dm[i]), L.mul(Fzi[i], dm[j]), s=-1)
                t2 = L.add(L.mul(Fzi[j], Di[i]), L.mul(Fzi[i], Di[j]), s=-1)
                cols.append(L.add(L.mul(D, t1), L.shift(t2, m), s=-1))
    for m in L.monomials(4):
        cols.append(L.shift(Fz, m, -1))
    # rhs = -(D * div(N0) - N0(D)), N0 = A e_k: div = d_k A ; N0(D) = A * D_k
    rhs = L.add(L.mul(D, L.pderiv(Az, k)), L.mul(Az, Di[k]), s=-1)
    rhs = {e: -v for e, v in rhs.items()}
    rows = {}
    for P in cols + [rhs]:
        for e in P: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    M = fmpz_mat(R, C + 1)
    for j, P in enumerate(cols):
        for e, v in P.items(): M[rows[e], j] = v
    for e, v in rhs.items(): M[rows[e], C] = v
    ra = M.rank()
    M0 = fmpz_mat(R, C)
    for j, P in enumerate(cols):
        for e, v in P.items(): M0[rows[e], j] = v
    r0 = M0.rank()
    return {'D': f'dF/d{L.VARS[k]}', 'd': 5, 'rows': R, 'cols': C, 'rank': r0, 'rank_aug': ra, 'consistent': ra == r0, 'seconds': round(time.time() - t0, 1)}
for k in (0, 3):
    res = koszul_exact(k)
    L.log({'check': 'exact: Koszul-reduced E2, no parity restriction, integer model', **res,
           'result': 'inconsistent over Q => no certificate with this denominator (exact)' if not res['consistent'] else 'CONSISTENT?!'})
