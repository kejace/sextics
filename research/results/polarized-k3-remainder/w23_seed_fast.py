#!/usr/bin/env python3
"""Vectorized F_p enumeration of profile-B incidence seeds on W_{23,f} (see w23_seed_search.py for conventions).
x = kx t ; a = c0 + c1 t ; y = k0 + c3 t ; S = k1 + S1 t + ... + S6 t^6 ; T = k2 + T1 t + ... + T4 t^4.
H = 7 f(x,a,y) + 7 S^2 + 7 T^3, 7f = 7x^6 + 5x^5 a + 5x^2 a^3 y + 15 x y^5 + 7a^6 + 10 a^5 y + 7 y^6.
Degrees 7..12 involve only (S1..S6, T1..T4, k2): solved triangularly; degree 6 fixes k1; degrees 0..5 are checked
over all (c0,c1,k0,c3) using a precomputed table of f-coefficients.  Then Jacobian / freeness matrices."""
import sys, json, time, itertools, argparse
import numpy as np
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from w23_seed_search import (F7, DF, evalpoly_coeffs, polymul, polypow, padto, grad_on_curve, jacobian,
                             freeness_matrix, det_mod_p, rank_mod_p)

def run(p, kx, out, maxseeds=10**9):
    t0 = time.time()
    X = np.array([0, kx % p], dtype=np.int64)
    lows = list(itertools.product(range(p), repeat=4))        # (c0,c1,k0,c3)
    ftab = np.zeros((len(lows), 13), dtype=np.int64)
    for i,(c0,c1,k0,c3) in enumerate(lows):
        ftab[i] = evalpoly_coeffs(F7, X, np.array([c0,c1]), np.array([k0,c3]), p)
    slots = [('a',0),('a',1),('y',1)] + [('S',k) for k in range(1,7)] + [('T',k) for k in range(1,5)]
    found = 0; good = 0; tops = 0
    fo = open(out, 'w')
    for T1,T2,T3,T4,k2 in itertools.product(range(p), repeat=5):
        if T4 == 0: continue
        T = np.array([k2,T1,T2,T3,T4], dtype=np.int64)
        t3 = padto(polypow(T,3,p), 13)
        rhs = (-T4**3) % p
        for S6 in [s for s in range(1,p) if (s*s) % p == rhs]:
            S = np.zeros(7, dtype=np.int64); S[6] = S6
            inv14S6 = pow((14*S6) % p, p-2, p)
            for deg in range(11, 6, -1):
                s2 = padto(polymul(S,S,p),13)
                cur = (7*s2[deg] + 7*t3[deg]) % p
                S[deg-6] = (-cur*inv14S6) % p
            tops += 1
            s2_part6 = sum(int(S[i])*int(S[6-i]) for i in range(1,6)) % p
            k1 = (-(ftab[:,6] + 7*s2_part6 + 7*t3[6]) * inv14S6) % p
            ok = np.ones(len(lows), dtype=bool)
            for j in range(0,6):
                inner = sum(int(S[i])*int(S[j-i]) for i in range(1,j)) % p
                s2j = (2*k1*int(S[j]) + inner) % p if j >= 1 else (k1*k1) % p
                val = (ftab[:,j] + 7*s2j + 7*t3[j]) % p
                ok &= (val == 0)
            for idx in np.nonzero(ok)[0]:
                c0,c1,k0,c3 = lows[idx]
                Sfull = S.copy(); Sfull[0] = int(k1[idx])
                A_ = np.array([c0,c1]); Y_ = np.array([k0,c3])
                g = grad_on_curve(X, A_, Y_, Sfull, T, p)
                J = jacobian(g, slots, p); dJ = det_mod_p(J, p)
                rF = rank_mod_p(freeness_matrix(g, p, True), p); rF2 = rank_mod_p(freeness_matrix(g, p, False), p)
                rec = {'p':p,'profile':'B','kx':kx,'k':[int(k0),int(k1[idx]),int(k2)],'a':[int(c0),int(c1)],'y':[int(k0),int(c3)],
                       'S':[int(v) for v in Sfull],'T':[int(v) for v in T],'detJ':int(dJ),'rank_free16':int(rF),'rank_free14_nox':int(rF2),
                       'gradS0':int(g['S'][0]),'gradT0':int(g['T'][0]),'grady0':int(g['y'][0]),'grada0':int(g['a'][0])}
                fo.write(json.dumps(rec)+'\n'); found += 1
                if dJ and rF == 12: good += 1
                if found >= maxseeds: break
    fo.close()
    return {'p':p,'kx':kx,'top_solutions':tops,'seeds':found,'seeds_detJ_nonzero_and_free':good,'seconds':round(time.time()-t0,1)}

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--p',type=int,default=5); ap.add_argument('--kx',type=int,default=4); ap.add_argument('--out',required=True)
    a = ap.parse_args(); print(json.dumps(run(a.p,a.kx,a.out)))
