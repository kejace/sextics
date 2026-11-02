#!/usr/bin/env python3
"""Random vectorized search for profile-A incidence seeds on W_{23,f} over F_p for larger p (11, 13, ...).
Same conventions as w23_seed_fastA.py (x = kx t; a,y quadratic; S deg 6 with S6 = param; T deg 4 with T4 = param).
For each random batch of (a0,a1,a2,k0,y1,y2,T0..T4) we solve S6 (sqrt), S5..S0 triangularly from degrees 12..6 and test
degrees 0..5.  Successful seeds are checked for invertible incidence Jacobian and freeness rank and written as JSON lines."""
import sys, json, time, argparse
import numpy as np
sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from w23_seed_search import (F7, grad_on_curve, jacobian, freeness_matrix, det_mod_p, rank_mod_p, polypow, padto)

def batch_f_coeffs(X1, A, Y, p):
    """A: (n,3), Y: (n,3) integer arrays; x = X1 t. Return (n,13) coefficient array of 7f(x,a,y) mod p (vectorized)."""
    n = A.shape[0]
    def pw(P, k):   # P: (n, d+1) coefficient array -> P^k
        R = np.zeros((n,1), dtype=np.int64); R[:,0] = 1
        for _ in range(k):
            R2 = np.zeros((n, R.shape[1]+P.shape[1]-1), dtype=np.int64)
            for i in range(R.shape[1]):
                for j in range(P.shape[1]):
                    R2[:,i+j] = (R2[:,i+j] + R[:,i]*P[:,j]) % p
            R = R2
        return R
    Xp = np.zeros((n,2), dtype=np.int64); Xp[:,1] = X1
    out = np.zeros((n,13), dtype=np.int64)
    cache = {}
    def P(name, k):
        if (name,k) not in cache:
            base = {'x':Xp,'a':A,'y':Y}[name]
            cache[(name,k)] = pw(base, k)
        return cache[(name,k)]
    for (i,j,k), c in F7.items():
        term = np.zeros((n,1), dtype=np.int64); term[:,0] = c % p
        for name, e in (('x',i),('a',j),('y',k)):
            Q = P(name, e)
            T2 = np.zeros((n, term.shape[1]+Q.shape[1]-1), dtype=np.int64)
            for u in range(term.shape[1]):
                for v in range(Q.shape[1]):
                    T2[:,u+v] = (T2[:,u+v] + term[:,u]*Q[:,v]) % p
            term = T2
        out[:, :term.shape[1]] = (out[:, :term.shape[1]] + term) % p
    return out

def run(p, kx, out, trials, batch=200000, seed=1):
    t0 = time.time(); rng = np.random.default_rng(seed)
    sq = {}
    for s in range(1,p): sq.setdefault((s*s)%p, []).append(s)
    sqrt1 = np.zeros(p, dtype=np.int64)   # one sqrt or 0
    for v,l in sq.items(): sqrt1[v] = l[0]
    inv7 = pow(7 % p, p-2, p)
    slots = [('a',0),('a',1),('a',2),('y',1),('y',2)] + [('S',k) for k in range(1,6)] + [('T',k) for k in range(1,4)]
    found = good = 0; done = 0
    fo = open(out, 'a')
    while done < trials:
        n = min(batch, trials - done); done += n
        A = rng.integers(0, p, size=(n,3)); Y = rng.integers(0, p, size=(n,3))
        T = rng.integers(0, p, size=(n,5)); T[:,4] = rng.integers(1, p, size=n)
        f = batch_f_coeffs(np.full(n, kx % p), A, Y, p)
        # T^3 coefficients
        t3 = np.zeros((n,13), dtype=np.int64)
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    t3[:,i+j+k] = (t3[:,i+j+k] + T[:,i]*T[:,j]*T[:,k]) % p
        r = (-(inv7*f[:,12] + t3[:,12])) % p
        S6 = sqrt1[r]
        valid = S6 != 0
        for sign in (0,1):
            S6s = S6 if sign == 0 else (-S6) % p
            Sm = np.zeros((n,7), dtype=np.int64); Sm[:,6] = S6s
            inv14 = np.array([pow(int(14*s) % p, p-2, p) if s else 0 for s in S6s], dtype=np.int64)
            for deg in range(11, 5, -1):
                idx = deg-6; Sm[:,idx] = 0
                inner = np.zeros(n, dtype=np.int64)
                for i in range(7):
                    j = deg - i
                    if 0 <= j <= 6: inner = (inner + Sm[:,i]*Sm[:,j]) % p
                cur = (f[:,deg] + 7*inner + 7*t3[:,deg]) % p
                Sm[:,idx] = (-(cur)*inv14) % p
            Ssq = np.zeros((n,13), dtype=np.int64)
            for i in range(7):
                for j in range(7):
                    Ssq[:,i+j] = (Ssq[:,i+j] + Sm[:,i]*Sm[:,j]) % p
            H = (f + 7*Ssq + 7*t3) % p
            ok = valid & np.all(H == 0, axis=1)
            for idx_row in np.nonzero(ok)[0]:
                a_ = A[idx_row]; y_ = Y[idx_row]; S = Sm[idx_row].copy(); Tv = T[idx_row]
                X = np.array([0, kx % p], dtype=np.int64)
                g = grad_on_curve(X, a_, y_, S, Tv, p)
                J = jacobian(g, slots, p); dJ = det_mod_p(J, p)
                rF = rank_mod_p(freeness_matrix(g, p, True), p); rF2 = rank_mod_p(freeness_matrix(g, p, False), p)
                rec = {'p':p,'profile':'A','kx':kx,'params':{'k0':int(y_[0]),'k1':int(S[0]),'k2':int(Tv[0]),'S6':int(S[6]),'T4':int(Tv[4])},
                       'a':[int(v) for v in a_],'y':[int(v) for v in y_],'S':[int(v) for v in S],'T':[int(v) for v in Tv],
                       'detJ':int(dJ),'rank_free16':int(rF),'rank_free14_nox':int(rF2)}
                fo.write(json.dumps(rec)+'\n'); fo.flush(); found += 1
                if dJ and rF == 12: good += 1
        if (done // batch) % 25 == 0:
            print(json.dumps({'progress':done,'found':found,'good':good,'seconds':round(time.time()-t0,1)}), flush=True)
    fo.close()
    return {'p':p,'kx':kx,'trials':trials,'seeds':found,'good':good,'seconds':round(time.time()-t0,1)}

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--p',type=int,default=11); ap.add_argument('--kx',type=int,default=1)
    ap.add_argument('--out',required=True); ap.add_argument('--trials',type=int,default=10**6); ap.add_argument('--seed',type=int,default=1)
    a = ap.parse_args(); print(json.dumps(run(a.p,a.kx,a.out,a.trials,seed=a.seed)))
