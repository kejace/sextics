"""Served-target spaces (Y-odd, q=1) for named denominators of one lane, mod p, compressed rows.
Usage: python3 tg.py <lane> <batchname> [margin]
Candidates are defined in cands(lane) below.  Appends JSONL to /tmp/hodge/results/_theory-generic/served.jsonl
and to the lane log.  Each record: lane, name, d, parity(D), served_dim (incl. the 18 Jacobian classes), served_mod_J,
A_served (is the lane target in the image), seconds.
"""
import sys, time, json, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d, jacobian
from hodge.core import mul, add, sub, scale, shift, monomials, power, parse
import numpy as np
from flint import nmod_mat
from fractions import Fraction

WD = '/tmp/hodge/results/_theory-generic'
OUT = os.path.join(WD, 'served.jsonl')

def served_space(F, q, D, parity=1, p=None, margin=60, seed=3):
    p = p or core.PRIMES[0]
    A0 = {monomials(6 * q, parity)[0]: 1}
    S = fixed_d.FixedDSystem(F, A0, q, D, p, with_target_space=False, verbose=False)
    Dp = core.reduce_mod(D, p); D2 = mul(Dp, Dp, p)
    amons = monomials(6 * q, parity)
    for m in amons:
        S.cols.append(('A', -1, m)); S.colpolys.append(shift(D2, m, 1, p))
    for col in S.colpolys[S.ncols:]:
        for e in col:
            if e not in S.rows:
                S.rows[e] = len(S.rows)
    S.nrows = len(S.rows); S.ncols = len(S.cols)
    M = fixed_d._sparse_matrix(S, with_rhs=False)
    C = S.ncols
    k = min(S.nrows, C + margin)
    Cm = fixed_d._compress(M, p, k, seed)
    Mn = nmod_mat(k, C, Cm.ravel().tolist(), p)
    del Cm
    X, nullity = Mn.nullspace()
    a_idx = list(range(C - len(amons), C))
    if nullity == 0:
        return 0, [], S.nrows, C
    proj = nmod_mat(len(a_idx), nullity, [int(X[j, kk]) for j in a_idx for kk in range(nullity)], p)
    T = proj.transpose(); Tr, rk = T.rref()
    basis = []
    for i in range(rk):
        poly = {}
        for jj, m in enumerate(amons):
            v = int(Tr[i, jj])
            if v: poly[m] = v
        basis.append(poly)
    return rk, basis, S.nrows, C

def rank_polys(polys, mons, p):
    col = {e: j for j, e in enumerate(mons)}
    flat = [0] * (len(polys) * len(mons))
    for r, P in enumerate(polys):
        for e, v in P.items():
            flat[r * len(mons) + col[e]] = v % p
    return nmod_mat(len(polys), len(mons), flat, p).rank()

def analyse(F, A, q, D, p=None, margin=60):
    p = p or core.PRIMES[0]
    rk, basis, rows, cols = served_space(F, q, D, 1, p, margin)
    mons = monomials(6, 1)
    J = []
    Fp = core.reduce_mod(F, p)
    for i in range(6):
        Fi = core.pderiv(Fp, i, p)
        for j in range(6):
            e = [0]*6; e[j] = 1
            g = shift(Fi, tuple(e), 1, p)
            if core.parity_of(g) == 1: J.append(g)
    rJ = rank_polys(J, mons, p)
    rJB = rank_polys(J + basis, mons, p)
    rJA = rank_polys(J + basis + [core.reduce_mod(A, p)], mons, p)
    return dict(served_dim=rk, J_rank=rJ, served_mod_J=rJB - rJ, A_served=(rJA == rJB), rows=rows, cols=cols)

X0, X1, X2, Y0, Y1, Y2 = [parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
def pw(P, k): return power(P, k)
def ssum(*ps):
    r = {}
    for P in ps: r = add(r, P)
    return r
def pderivs(F): return [core.pderiv(F, i) for i in range(6)]
def E_X(F): Fi = pderivs(F); return ssum(mul(X0, Fi[0]), mul(X1, Fi[1]), mul(X2, Fi[2]))
def E_Y(F): Fi = pderivs(F); return ssum(mul(Y0, Fi[3]), mul(Y1, Fi[4]), mul(Y2, Fi[5]))
def bipart(F, kx):
    return {e: c for e, c in F.items() if e[0]+e[1]+e[2] == kx}

def cands(lane, F):
    C = []
    if lane == 'short-power-presentation':
        L = parse('X0 + 2*X1 + 3*X2'); M = parse('Y0 + 3*Y1 + 5*Y2')
        eX = parse('X0**6 + X1**6 + X2**6'); eY = parse('Y0**6 + Y1**6 + Y2**6')
        c = Fraction(5, 7**7)
        U6 = add(pw(add(L, M), 6), pw(sub(L, M), 6))   # (L+M)^6 + (L-M)^6, even
        V6 = sub(pw(add(L, M), 6), pw(sub(L, M), 6))   # odd
        C += [('vertexline_X0X1Y0Y1', 'D_I |I|=4: cone vertex line {X0=X1=Y0=Y1=0} over Fermat-type sextic surface', ssum(pw(X0,6),pw(X1,6),pw(Y0,6),pw(Y1,6))),
              ('eX', 'cone vertex plane {X=0} over Fermat curve', eX),
              ('eX_plus_2cL6', 'F restricted to M=0 (cone vertex line {X=0,L=0}? no: depends on X only) ', add(eX, scale(pw(L,6), 2*c))),
              ('U6_coupling', '(L+M)^6+(L-M)^6 = union of 6 hyperplane sections through P^3 {L=M=0}', U6),
              ('V6_coupling_odd', '(L+M)^6-(L-M)^6 odd, 6 hyperplanes through {L=M=0}', V6),
              ('X0^6+Y0^6', 'binary block, 6 hyperplanes through P^3 {X0=Y0=0}', add(pw(X0,6), pw(Y0,6))),
              ('X0^6+X1^6+Y0^6', '|I|=3 mixed: cone vertex plane {X0=X1=Y0=0}', ssum(pw(X0,6),pw(X1,6),pw(Y0,6))),
              ('EX_F', 'polar: X-Euler derivative, tangency locus', E_X(F)),
              ]
    if lane == 'eight-plane-configuration':
        Q = [add(pw(X0,2),pw(Y0,2)), add(pw(X1,2),pw(Y1,2)), add(pw(X2,2),pw(Y2,2))]
        Fi = pderivs(F)
        R0 = sub(mul(Y0, Fi[0]), mul(X0, Fi[3]))  # odd, block-0 rotation field applied to F; vanishes on all 8 planes
        C += [('R0_F_odd', 'rotation field Y0 d/dX0 - X0 d/dY0 applied to F (odd, d=6); tangency locus contains the 8 planes', R0),
              ('Q0Q1Q2', 'reduced equation of the 8-plane arrangement; union of 6 hyperplane sections', mul(mul(Q[0],Q[1]),Q[2])),
              ('Q0^2Q1', 'double along 4 planes', mul(pw(Q[0],2), Q[1])),
              ('Q0^3', 'triple along V(Q0)', pw(Q[0],3)),
              ('vertexline_X0X1Y0Y1', 'cone vertex line over Fermat-type sextic surface', ssum(pw(X0,6),pw(X1,6),pw(Y0,6),pw(Y1,6))),
              ('EX_F', 'polar X-Euler', E_X(F)),
              ]
    if lane == 'quadratic-cubic-pullback':
        QS = ["X0**2 + 5*X1*X2/7 + 5*Y0*Y1/7", "5*X0*X2/11 + X1**2 + 5*Y1*Y2/11", "5*X0*X1/13 + X2**2 + 5*Y0*Y2/13",
              "5*X0*X1/17 + Y0**2 + 5*Y1*Y2/17", "5*X1*X2/19 + 5*Y0*Y2/19 + Y1**2", "5*X0*X2/23 + 5*Y0*Y1/23 + Y2**2"]
        Q = [parse(s) for s in QS]
        C += [('Q0Q1Q2', 'pullback of 3 hyperplane sections of the cubic', mul(mul(Q[0],Q[1]),Q[2])),
              ('Q0^3+Q1^3', 'pullback of the Eckardt-type reducible cubic section (3 hyperplanes through a P^3)', add(pw(Q[0],3), pw(Q[1],3))),
              ('Q0^3+Q1^3+Q2^3+Q3^3', 'pullback of a 4-term Fermat cubic section (cone over cubic surface... vertex line)', ssum(pw(Q[0],3),pw(Q[1],3),pw(Q[2],3),pw(Q[3],3))),
              ('vertexline_X0X1Y0Y1', 'cone vertex line over Fermat-type sextic surface', ssum(pw(X0,6),pw(X1,6),pw(Y0,6),pw(Y1,6))),
              ('F60', 'F_{6,0}: cone vertex plane {X=0}', bipart(F, 6)),
              ('F60+F06', 'F_{6,0}+F_{0,6}', add(bipart(F,6), bipart(F,0))),
              ]
    if lane == 'dense-rational':
        C += [('F60', 'F_{6,0}=a(X): cone vertex plane P_Y over the plane sextic a', bipart(F, 6)),
              ('F60+F06', 'a(X)+d(Y): join-type', add(bipart(F,6), bipart(F,0))),
              ('vertexline_X0X1Y0Y1', 'cone vertex line over Fermat-type sextic surface', ssum(pw(X0,6),pw(X1,6),pw(Y0,6),pw(Y1,6))),
              ('X0^6+Y0^6', 'binary block: 6 hyperplanes through P^3', add(pw(X0,6), pw(Y0,6))),
              ('EY_F', 'polar Y-Euler', E_Y(F)),
              ]
    return C

if __name__ == '__main__':
    lane = sys.argv[1]; batch = sys.argv[2]
    margin = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    only = sys.argv[4].split(',') if len(sys.argv) > 4 else None
    F, A, q, rec = core.load_instance(lane)
    for name, note, D in cands(lane, F):
        if only and name not in only: continue
        d = core.degree(D); par = core.parity_of(D)
        t0 = time.time()
        try:
            ok, info = fixed_d.e1_solvable(F, A, q, D)
            r = analyse(F, A, q, D, margin=margin)
            r.update(dict(lane=lane, batch=batch, name=name, note=note, d=d, parity_D=par, e1_for_A=bool(ok), seconds=round(time.time()-t0,1), engine=2, ts=time.strftime('%Y-%m-%dT%H:%M:%S')))
        except Exception as ex:
            r = dict(lane=lane, batch=batch, name=name, d=d, parity_D=par, error=repr(ex), seconds=round(time.time()-t0,1))
        print(json.dumps(r), flush=True)
        with open(OUT, 'a') as f: f.write(json.dumps(r) + '\n')
        with open(f'/tmp/hodge/results/{lane}/log.jsonl', 'a') as f:
            f.write(json.dumps(dict(D=name, d=d, family='theory-generic-served', served_dim=r.get('served_dim'), served_mod_J=r.get('served_mod_J'), A_served=r.get('A_served'), e1=r.get('e1_for_A'), seconds=r.get('seconds'), note=note, engine=2, agent='theory-generic')) + '\n')
