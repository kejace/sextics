#!/usr/bin/env python3
"""Profile A (dimension-matched, 5 parameters) incidence seeds on W_{23,f} over F_p:
x = kx t ; a = a0 + a1 t + a2 t^2 ; y = k0 + y1 t + y2 t^2 ; S = k1 + S1 t + ... + S5 t^5 + k3 t^6 ; T = k2 + T1 t + T2 t^2 + T3 t^3 + k4 t^4.
Parameters (k0,k1,k2,k3,k4) = constant terms of y,S,T and the top coefficients S6,T4.  Unknowns (13): a0,a1,a2,y1,y2,S1..S5,T1,T2,T3.
Reason for the extra two parameters: on W_{23} (-K = O(2)) weighted conics through a point form a 2-dim family, so the
3-parameter normalization of the crossed-ternary construction leaves a positive-dimensional incidence (det J = 0, see profile B).
H = 7 f + 7 S^2 + 7 T^3.  Degrees 12..6 solved triangularly for S6(root),S5,...,S1,S0=k1; degrees 0..5 checked (vectorized)."""
import sys, json, time, itertools, argparse
import numpy as np
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from w23_seed_search import (F7, evalpoly_coeffs, polymul, polypow, padto, grad_on_curve, jacobian, freeness_matrix, det_mod_p, rank_mod_p)

def sqrt_table(p):
    tab = {}
    for s in range(1, p): tab.setdefault((s*s) % p, []).append(s)
    return tab

def run(p, kx, out, maxseeds=10**9, want_good=50):
    t0 = time.time()
    X = np.array([0, kx % p], dtype=np.int64)
    ay = list(itertools.product(range(p), repeat=6))   # a0,a1,a2,k0,y1,y2
    n = len(ay)
    ftab = np.zeros((n, 13), dtype=np.int64)
    for i,(a0,a1,a2,k0,y1,y2) in enumerate(ay):
        ftab[i] = evalpoly_coeffs(F7, X, np.array([a0,a1,a2]), np.array([k0,y1,y2]), p)
    sq = sqrt_table(p)
    slots = [('a',0),('a',1),('a',2),('y',1),('y',2)] + [('S',k) for k in range(1,6)] + [('T',k) for k in range(1,4)]
    found = good = 0; tops = 0
    fo = open(out, 'w')
    inv7 = pow(7 % p, p-2, p)
    for T1,T2,T3,T4,k2 in itertools.product(range(p), repeat=5):
        if T4 == 0: continue
        T = np.array([k2,T1,T2,T3,T4], dtype=np.int64)
        t3 = padto(polypow(T,3,p), 13)
        tops += 1
        # degree 12: 7 f12 + 7 S6^2 + 7 t3_12 = 0 -> S6^2 = -(f12 + t3_12)
        r = (-(inv7*ftab[:,12] + t3[12])) % p
        for sign in (0,1):
            S6 = np.array([ (sq[v][sign] if (v in sq and len(sq[v])>sign) else 0) for v in r ], dtype=np.int64)
            valid = S6 != 0
            if not valid.any(): continue
            Sm = np.zeros((n,7), dtype=np.int64); Sm[:,6] = S6
            inv14S6 = np.array([pow(int(14*s) % p, p-2, p) if s else 0 for s in S6], dtype=np.int64)
            # degrees 11..6: [t^deg]: 7( f_deg + (S^2)_deg + t3_deg ) = 0 ; (S^2)_deg = 2 S6 S_{deg-6} + sum_{i+j=deg, i,j<6, i>=j... } ...
            for deg in range(11, 5, -1):
                idx = deg - 6
                Sm[:,idx] = 0
                inner = np.zeros(n, dtype=np.int64)
                for i in range(0, 7):
                    j = deg - i
                    if 0 <= j <= 6:
                        inner = (inner + Sm[:,i]*Sm[:,j]) % p
                # equation: ftab[deg] + 7*(2 S6 S_idx + inner) + 7 t3[deg] = 0  ->  S_idx = -(ftab[deg] + 7 inner + 7 t3[deg]) * inv(14 S6)
                cur = (ftab[:,deg] + 7*inner + 7*t3[deg]) % p
                Sm[:,idx] = (-(cur) * inv14S6) % p
            # now S0..S6 set (S0 = k1). check degrees 0..5
            Ssq = np.zeros((n,13), dtype=np.int64)
            for i in range(7):
                for j in range(7):
                    Ssq[:,i+j] = (Ssq[:,i+j] + Sm[:,i]*Sm[:,j]) % p
            H = (ftab + 7*Ssq + 7*t3[None,:]) % p
            ok = valid & np.all(H == 0, axis=1)
            for idx_row in np.nonzero(ok)[0]:
                a0,a1,a2,k0,y1,y2 = ay[idx_row]
                S = Sm[idx_row].copy()
                A_ = np.array([a0,a1,a2]); Y_ = np.array([k0,y1,y2])
                g = grad_on_curve(X, A_, Y_, S, T, p)
                J = jacobian(g, slots, p); dJ = det_mod_p(J, p)
                rF = rank_mod_p(freeness_matrix(g, p, True), p); rF2 = rank_mod_p(freeness_matrix(g, p, False), p)
                rec = {'p':p,'profile':'A','kx':kx,'params':{'k0':int(k0),'k1':int(S[0]),'k2':int(k2),'S6':int(S[6]),'T4':int(T4)},
                       'a':[int(a0),int(a1),int(a2)],'y':[int(k0),int(y1),int(y2)],'S':[int(v) for v in S],'T':[int(v) for v in T],
                       'detJ':int(dJ),'rank_free16':int(rF),'rank_free14_nox':int(rF2),
                       'gradS0':int(g['S'][0]),'gradT0':int(g['T'][0]),'grady0':int(g['y'][0]),'grada0':int(g['a'][0])}
                fo.write(json.dumps(rec)+'\n'); found += 1
                if dJ and rF == 12: good += 1
        if found >= maxseeds: break
    fo.close()
    return {'p':p,'kx':kx,'T_choices':tops,'seeds':found,'seeds_detJ_nonzero_and_free':good,'seconds':round(time.time()-t0,1)}

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--p',type=int,default=5); ap.add_argument('--kx',type=int,default=4); ap.add_argument('--out',required=True)
    a = ap.parse_args(); print(json.dumps(run(a.p,a.kx,a.out)))
