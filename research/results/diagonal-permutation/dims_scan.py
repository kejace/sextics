"""Dimension scan of J_d^+ ∩ (linear ideal)^k for S3/tau-adapted subspaces, d = 7,8,9 (mod p, fast)."""
import sys, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from jideal import find
from fractions import Fraction as Fr
I3 = [[1,0,0],[0,1,0],[0,0,1]]
H = Fr(1,2)
# sigma-adapted: X = Mx z with z=(p,r,a): X0=(p+r)/2, X1=(p-r)/2, X2=a ; Y similarly (q,s,b)
Msig = [[H,H,0],[H,-H,0],[0,0,1]]
# 3-cycle adapted: z = (e, r1, r2) with X0 = e + r1, X1 = e + r2, X2 = e - r1 - r2  (so r1 = (2X0 - X1 - X2)/3 etc.)
Mcyc = [[1,1,0],[1,0,1],[1,-1,-1]]
cases = [
 ('orig: (X0,Y0,X1)^3 [Fermat-type cone plane for T0]', I3, I3, [0,3,1], 3),
 ('orig: (X0,Y0,X1)^4', I3, I3, [0,3,1], 4),
 ('orig: (X0,Y0,Y1)^3', I3, I3, [0,3,4], 3),
 ('orig: (X0,Y0,Y1)^4', I3, I3, [0,3,4], 4),
 ('orig: (X0,Y0,X2)^3', I3, I3, [0,3,2], 3),
 ('orig: (X0,Y0,Y2)^3', I3, I3, [0,3,5], 3),
 ('orig: (X2,Y2)^2 [non-normal along X∩{X2=Y2=0}]', I3, I3, [2,5], 2),
 ('orig: (X0,Y0)^2', I3, I3, [0,3], 2),
 ('orig: (X0,X1)^2', I3, I3, [0,1], 2),
 ('orig: (Y0,Y1)^2', I3, I3, [3,4], 2),
 ('orig: (X0,Y0)^3', I3, I3, [0,3], 3),
 ('orig: (X2,Y2)^3', I3, I3, [2,5], 3),
 ('orig: (X0,X1,X2,Y0)^4 [vertex line {X=0,Y0=0}]', I3, I3, [0,1,2,3], 4),
 ('orig: (X0,X1,Y0,Y1)^4 [vertex line = (X2,Y2)-line]', I3, I3, [0,1,3,4], 4),
 ('sigma: (r,s)^2 [non-normal along the (01)-fixed surface]', Msig, Msig, [1,4], 2),
 ('sigma: (r,s)^3', Msig, Msig, [1,4], 3),
 ('sigma: (p,q,a,b)^4 [quadruple points on the anti-invariant line]', Msig, Msig, [0,3,2,5], 4),
 ('sigma: (p,q,a,b)^3', Msig, Msig, [0,3,2,5], 3),
 ('sigma: (r,s,a)^3 [cone plane {X0=X1,Y0=Y1,X2=0}]', Msig, Msig, [1,4,2], 3),
 ('sigma: (r,s,b)^3', Msig, Msig, [1,4,5], 3),
 ('sigma: (p,q,a)^3', Msig, Msig, [0,3,2], 3),
 ('sigma: (p,q,b)^3', Msig, Msig, [0,3,5], 3),
 ('sigma: (r,a,s,b)^4 [vertex line {X0=X1,X2=0,Y0=Y1,Y2=0}... = fixed line of (01)? no: p,q free]', Msig, Msig, [1,2,4,5], 4),
 ('cyc: (r1,r2,s1,s2)^4 [quadruple at the 6 S3-fixed points]', Mcyc, Mcyc, [1,2,4,5], 4),
 ('cyc: (r1,r2,s1,s2)^3', Mcyc, Mcyc, [1,2,4,5], 3),
 ('cyc: (r1,r2)^2 [non-normal along X∩{X0=X1=X2}]', Mcyc, Mcyc, [1,2], 2),
 ('cyc: (s1,s2)^2 [non-normal along X∩{Y0=Y1=Y2}]', Mcyc, Mcyc, [4,5], 2),
 ('cyc: (r1,r2,e)^3 = (X0,X1,X2)^3 check', Mcyc, Mcyc, [0,1,2], 3),
]
for d in (7, 8, 9):
    for name, Mx, My, idx, k in cases:
        t0 = time.time()
        try:
            r = find(Mx, My, idx, k, d, verbose=False)
            print(f'd={d} {name}: nontrivial dim = {r["nontrivial"]}  ({time.time()-t0:.1f}s)', flush=True)
        except Exception as ex:
            print(f'd={d} {name}: ERROR {ex}', flush=True)
