"""Replay an exact Farkas certificate: python3 verify_farkas.py farkas_<name>.pkl
Rebuilds nothing from the folded system: it only uses (F, A, D, y) and checks y . column = 0 for every unknown of the
full (unreduced) fixed-D system and y . rhs != 0, in exact integer arithmetic."""
import sys, pickle, time
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import numpy as np
from fractions import Fraction
from hodge import core
import charsolve as CS
rec = pickle.load(open(sys.argv[1], 'rb'))
F, A, q, _ = core.load_instance('symmetric-sextic-residual')
lf = rec['scaleF']; Fz = {e: int(Fraction(c) * lf) for e, c in F.items()}
Az = {e: int(Fraction(c)) for e, c in A.items()}
D = {tuple(int(x) for x in k.strip('()').split(',')): Fraction(v) for k, v in rec['D'].items()}
ld = rec['scaleD']; Dz = {e: int(c * ld) for e, c in D.items()}
d = core.degree(Dz); pD = core.parity_of(Dz); pA = 1
degN, degQ, degE = d + 1, 2 * d - 6, 2 * d + 6
G = CS.Group([[0, 1, 2], [3, 4, 5]], {'triv': [1, 1], 'sgnX': [-1, 1], 'sgnY': [1, -1], 'sgnsgn': [-1, -1]}[rec['chi']])
codes = np.array(rec['orbit_codes'], dtype=np.int64); y = rec['y_orbit']
rowmons = core.monomials(degE, pA)
idx = np.searchsorted(codes, G.canon(np.array(rowmons, dtype=np.int64)))
yfull = {m: y[i] for m, i in zip(rowmons, idx.tolist()) if y[i]}
def dot(P): return sum(yfull.get(e, 0) * int(c) for e, c in P.items())
Fi = [core.pderiv(Fz, i) for i in range(6)]; Di = [core.pderiv(Dz, i) for i in range(6)]; FD = core.mul(Fz, Dz)
t0 = time.time(); bad = 0; ncols = 0
for s in range(6):
    G_s = core.sub(core.scale(core.mul(Dz, Fi[s]), -1), core.mul(Fz, Di[s]))   # coefficient of N_s-monomial m:  (-D F_s - F D_s) m + F D d(m)/dz_s
    for m in core.monomials(degN, CS.pN_of(pA, pD, s)):
        P = core.shift(G_s, m, 1)
        if m[s]:
            m2 = list(m); m2[s] -= 1
            P = core.add(P, core.shift(FD, tuple(m2), m[s]))
        if dot(P) != 0: bad += 1
        ncols += 1
if degQ >= 0:
    F2 = core.scale(core.mul(Fz, Fz), -1)
    for m in core.monomials(degQ, pA):
        if dot(core.shift(F2, m, 1)) != 0: bad += 1
        ncols += 1
yb = dot(core.scale(core.mul(Az, core.mul(Dz, Dz)), -1))
print(f"{rec['label']}: d={d} full-system columns checked={ncols}, violations={bad}, y.rhs={yb} ({time.time()-t0:.1f}s)")
print('VERDICT:', 'inconsistent over Q (certificate valid)' if bad == 0 and yb != 0 else 'CERTIFICATE INVALID')
