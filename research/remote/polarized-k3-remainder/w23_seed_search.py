#!/usr/bin/env python3
"""Wave-2 construction lens: polynomial-curve incidence seeds on the second elliptic companion
   W_{23,f} = { H = 7 f(x,a,y) + 7 S^2 + 7 T^3 = 0 } in P(1,1,1,3,2), f the ternary sextic of polarized-k3-remainder
   (7f = 7x^6 + 5x^5 a + 5x^2 a^3 y + 15 x y^5 + 7a^6 + 10 a^5 y + 7 y^6, variables (x,a,y) = (z0,z1,z2)).
Curve of weighted degree e=2 with the distinguished coordinate x = k_x t (zero constant term), endpoint constants as
parameters, remaining coefficients as incidence unknowns; square 13x13 incidence system (degrees 0..12 of H(phi)).
Profiles:
  B: a = c0 + c1 t, y = k0 + c3 t (linear), S = k1 + S1 t + ... + S6 t^6, T = k2 + T1 t + ... + T4 t^4   (13 unknowns)
  A: a = c0 + c1 t + c2 t^2, y = k0 + c3 t + c4 t^2, S up to t^5 (S6 = k3 parameter), T up to t^3 (T4 = k4 parameter)
Brute force over F_p exploiting the triangular structure in the top degrees; checks incidence Jacobian (13x13) and the
freeness matrix (12 x 16; also the 12 x 14 version without the x-columns).  Output: JSON lines of seeds.
Usage: python3 w23_seed_search.py --p 5 --profile B --kx 4 --out seeds_B_p5.jsonl [--max 200]
"""
import argparse, itertools, json, sys, time
import numpy as np

def polymul(a, b, p):
    return np.convolve(a, b) % p

def polypow(a, n, p):
    r = np.array([1], dtype=np.int64)
    for _ in range(n):
        r = polymul(r, a, p)
    return r

def padto(a, n):
    out = np.zeros(n, dtype=np.int64); out[:len(a)] = a[:n] if len(a) >= n else a; return out

# 7f and its partials as dict monomial->coeff (integers)
F7 = {(6,0,0):7, (5,1,0):5, (2,3,1):5, (1,0,5):15, (0,6,0):7, (0,5,1):10, (0,0,6):7}
def evalpoly_coeffs(terms, X, A, Y, p, deg=12):
    """terms: dict (i,j,k)->c ; X,A,Y numpy coefficient arrays; return coefficient array of degree<=deg mod p."""
    out = np.zeros(deg+1, dtype=np.int64)
    pw = {}
    def pwr(v, n, key):
        if (key,n) not in pw: pw[(key,n)] = polypow(v, n, p)
        return pw[(key,n)]
    for (i,j,k), c in terms.items():
        t = np.array([c % p], dtype=np.int64)
        t = polymul(t, pwr(X,i,'x'), p); t = polymul(t, pwr(A,j,'a'), p); t = polymul(t, pwr(Y,k,'y'), p)
        out = (out + padto(t, deg+1)) % p
    return out

def partials():
    d = {'x':{}, 'a':{}, 'y':{}}
    for (i,j,k), c in F7.items():
        if i: d['x'][(i-1,j,k)] = d['x'].get((i-1,j,k),0) + c*i
        if j: d['a'][(i,j-1,k)] = d['a'].get((i,j-1,k),0) + c*j
        if k: d['y'][(i,j,k-1)] = d['y'].get((i,j,k-1),0) + c*k
    return d
DF = partials()

def H_coeffs(X, A, Y, S, T, p):
    f = evalpoly_coeffs(F7, X, A, Y, p)
    s2 = padto(polymul(S, S, p), 13); t3 = padto(polypow(T, 3, p), 13)
    return (f + 7*s2 + 7*t3) % p

def grad_on_curve(X, A, Y, S, T, p):
    g = {}
    g['x'] = evalpoly_coeffs(DF['x'], X, A, Y, p)
    g['a'] = evalpoly_coeffs(DF['a'], X, A, Y, p)
    g['y'] = evalpoly_coeffs(DF['y'], X, A, Y, p)
    g['S'] = padto((14*S) % p, 13)              # d/dS (7 S^2) = 14 S
    g['T'] = padto((21*polymul(T,T,p)) % p, 13)   # d/dT (7 T^3) = 21 T^2
    return g

def rank_mod_p(M, p):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i,c] % p: piv = i; break
        if piv is None: continue
        M[[r,piv]] = M[[piv,r]]
        inv = pow(int(M[r,c]), p-2, p)
        M[r] = (M[r]*inv) % p
        for i in range(rows):
            if i != r and M[i,c]: M[i] = (M[i] - M[i,c]*M[r]) % p
        r += 1
        if r == rows: break
    return r

def det_mod_p(M, p):
    M = M.copy() % p; n = M.shape[0]; det = 1
    for c in range(n):
        piv = None
        for i in range(c, n):
            if M[i,c] % p: piv = i; break
        if piv is None: return 0
        if piv != c: M[[c,piv]] = M[[piv,c]]; det = (-det) % p
        det = (det*int(M[c,c])) % p
        inv = pow(int(M[c,c]), p-2, p)
        for i in range(c+1, n):
            if M[i,c]: M[i] = (M[i] - M[i,c]*inv*M[c]) % p
    return det

def freeness_matrix(g, p, with_x=True):
    cols = []
    for coord, w in (('x',1),('a',1),('y',1),('S',3),('T',2)):
        if coord == 'x' and not with_x: continue
        for k in range(2*w):  # t^0..t^(w e -1), e=2
            v = np.zeros(12, dtype=np.int64)
            gc = g[coord]
            for i in range(12):
                if 0 <= i-k < len(gc): v[i] = gc[i-k] % p
            cols.append(v)
    return np.array(cols, dtype=np.int64).T  # 12 x ncols

def jacobian(g, slots, p):
    """slots: list of (coord, degree) for each unknown; rows = degrees 0..12 of H(phi)."""
    J = np.zeros((13, len(slots)), dtype=np.int64)
    for j,(coord,k) in enumerate(slots):
        gc = g[coord]
        for i in range(13):
            if 0 <= i-k < len(gc): J[i,j] = gc[i-k] % p
    return J

def search(p, profile, kx, maxseeds, out, log_every=200000):
    t0 = time.time(); found = 0; tried = 0
    X = np.array([0, kx % p], dtype=np.int64)
    fo = open(out, 'a')
    if profile == 'B':
        slots = [('a',0),('a',1),('y',1)] + [('S',k) for k in range(1,7)] + [('T',k) for k in range(1,5)]
        # enumerate (T1..T4, k2) ; solve S6..S1 from degrees 12..7 (S6 quadratic, rest linear: coefficient 14*S6)
        for T1,T2,T3,T4,k2 in itertools.product(range(p), repeat=5):
            if T4 == 0: continue
            T = np.array([k2,T1,T2,T3,T4], dtype=np.int64)
            t3 = padto(polypow(T,3,p), 13)
            # degree 12: 7 S6^2 + 7 T4^3 = 0 -> S6^2 = -T4^3
            rhs = (-T4**3) % p
            S6s = [s for s in range(1,p) if (s*s) % p == rhs]
            for S6 in S6s:
                S = np.zeros(7, dtype=np.int64); S[6] = S6
                ok = True
                for deg in range(11, 6, -1):   # solve S_{deg-6}
                    idx = deg - 6
                    # coefficient of t^deg in 7 S^2 + 7 T^3 (x,a,y part has degree<=6 in profile B)
                    s2 = padto(polymul(S,S,p),13)
                    cur = (7*s2[deg] + 7*t3[deg]) % p
                    # adding S_idx changes s2[deg] by 2*S_idx*S6 (idx+6 = deg) -> coefficient 14*S6
                    coef = (14*S6) % p
                    S[idx] = (-(cur) * pow(coef, p-2, p)) % p
                # now degrees 7..12 vanish; enumerate lower unknowns c0,c1,k0,c3,k1
                for c0,c1,k0,c3 in itertools.product(range(p), repeat=4):
                    A_ = np.array([c0,c1], dtype=np.int64); Y_ = np.array([k0,c3], dtype=np.int64)
                    f = evalpoly_coeffs(F7, X, A_, Y_, p)
                    # degree 6: 7(S^2)[6] + 7 t3[6] + f[6] = 0 with S0 = k1 unknown: (S^2)[6] = 2 k1 S6 + sum_{i=1..5} S_i S_{6-i}
                    s2_part = sum(int(S[i])*int(S[6-i]) for i in range(1,6)) % p
                    cur6 = (f[6] + 7*(s2_part) + 7*t3[6]) % p
                    k1 = (-(cur6) * pow((14*S6) % p, p-2, p)) % p
                    S[0] = k1
                    Hc = H_coeffs(X, A_, Y_, S, T, p)
                    tried += 1
                    if np.any(Hc % p):
                        S[0] = 0
                        continue
                    S[0] = k1
                    g = grad_on_curve(X, A_, Y_, S, T, p)
                    J = jacobian(g, slots, p); dJ = det_mod_p(J, p)
                    Fm = freeness_matrix(g, p, True); rF = rank_mod_p(Fm, p)
                    Fm2 = freeness_matrix(g, p, False); rF2 = rank_mod_p(Fm2, p)
                    rec = {'p':p,'profile':'B','kx':kx,'k':[int(k0),int(k1),int(k2)],
                           'a':[int(c0),int(c1)],'y':[int(k0),int(c3)],'S':[int(v) for v in S],'T':[int(v) for v in T],
                           'detJ':int(dJ),'rank_free16':int(rF),'rank_free14_nox':int(rF2),
                           'gradS_at_endpoint': int(g['S'][0]), 'gradT_at_endpoint': int(g['T'][0])}
                    fo.write(json.dumps(rec)+'\n'); fo.flush(); found += 1
                    S[0] = 0
                    if found >= maxseeds: fo.close(); return found, tried, time.time()-t0
                if tried % log_every < p**4: print(f'[{time.time()-t0:.0f}s] tried={tried} found={found}', flush=True)
    elif profile == 'A':
        slots = [('a',0),('a',1),('a',2),('y',1),('y',2)] + [('S',k) for k in range(1,6)] + [('T',k) for k in range(1,4)]
        # parameters k0 (y const), k1 (S const), k2 (T const), k3 = S6, k4 = T4 ; unknowns 13 ; brute force layered:
        # top degrees 12..: f(4t,a,y) now has degree 12 -> no triangular shortcut; use random sampling of the
        # 13 unknowns + 5 params is too big (5^18); instead enumerate (a2,y2,k3,k4) to satisfy degree-12 eq, then
        # sample the rest randomly.
        rng = np.random.default_rng(12345)
        samples = 0
        while found < maxseeds and samples < 4_000_000:
            vals = rng.integers(0, p, size=18); samples += 1; tried += 1
            c0,c1,c2,k0,c3,c4 = [int(v) for v in vals[:6]]
            Sv = [int(v) for v in vals[6:13]]   # S0..S6
            Tv = [int(v) for v in vals[13:18]]  # T0..T4
            A_ = np.array([c0,c1,c2]); Y_ = np.array([k0,c3,c4]); S = np.array(Sv); T = np.array(Tv)
            Hc = H_coeffs(X, A_, Y_, S, T, p)
            if np.any(Hc % p): continue
            g = grad_on_curve(X, A_, Y_, S, T, p)
            J = jacobian(g, slots, p); dJ = det_mod_p(J, p)
            Fm = freeness_matrix(g, p, True); rF = rank_mod_p(Fm, p)
            rec = {'p':p,'profile':'A','kx':kx,'params':{'k0':k0,'k1':Sv[0],'k2':Tv[0],'S6':Sv[6],'T4':Tv[4]},
                   'a':[c0,c1,c2],'y':[k0,c3,c4],'S':Sv,'T':Tv,'detJ':int(dJ),'rank_free16':int(rF)}
            fo.write(json.dumps(rec)+'\n'); fo.flush(); found += 1
    fo.close()
    return found, tried, time.time()-t0

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=int, default=5); ap.add_argument('--profile', default='B')
    ap.add_argument('--kx', type=int, default=4); ap.add_argument('--out', required=True)
    ap.add_argument('--max', type=int, default=100000)
    a = ap.parse_args()
    found, tried, secs = search(a.p, a.profile, a.kx, a.max, a.out)
    print(json.dumps({'done':True,'p':a.p,'profile':a.profile,'kx':a.kx,'found':found,'tried':tried,'seconds':round(secs,1)}))
