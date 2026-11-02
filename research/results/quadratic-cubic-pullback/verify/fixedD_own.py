"""Own-code fixed-denominator system (uncompressed, dense rank mod p) for (F, q=1, A).
E(N,Q) = A D^2 - q D sum N_i F_i + F (D sum d_i N_i - sum N_i d_i D) - F^(q+1) Q.
Functions: consistent(D, target) and served_dim(D). usage: python3 fixedD_own.py <case> [p]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback/verify')
import vpoly as V
from vlog import log
from flint import nmod_mat
p = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000007
F, A, q = V.load_instance_F_A()
chk = json.load(open('/home/kejace/devel/sextics/data/quadratic-cubic-pullback/checks.json'))
Q = [V.parse(s) for s in chk['quadrics']]
delta = V.parse(open('/tmp/hodge/results/quadratic-cubic-pullback/verify/delta_exact.txt').read())
F0 = V.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')

def build(Fx, D, target, with_target_cols=False):
    Fp, Dp = V.modp(Fx, p), V.modp(D, p)
    Tp = V.modp(target, p) if target else {}
    d = V.deg(D); pD = V.parity(D); pA = 1   # targets are Y-odd
    assert pD is not None
    Fi = [V.pderiv(Fp, i, p) for i in range(6)]
    Di = [V.pderiv(Dp, i, p) for i in range(6)]
    degN, degQ = d + 6*q - 5, 2*d - 6
    DF = [V.scale(V.mul(Dp, Fi[i], p), (-q) % p, p) for i in range(6)]     # -q D F_i
    FD = V.mul(Fp, Dp, p)
    FDi = [V.scale(V.mul(Fp, Di[i], p), p - 1, p) for i in range(6)]      # -F D_i
    Fq1 = V.scale(V.power(Fp, q + 1, p), p - 1, p)                        # -F^{q+1}
    cols, kinds = [], []
    for i in range(6):
        fpar = 0 if i < 3 else 1
        par_i = (pA - pD - fpar) % 2
        for m in V.monomials(degN, par_i):
            col = V.shift(DF[i], m, 1, p)
            col = V.add(col, V.shift(FDi[i], m, 1, p), p)
            if m[i]:
                m2 = list(m); m2[i] -= 1
                col = V.add(col, V.shift(FD, tuple(m2), m[i], p), p)       # F D d_i(m)
            cols.append(col); kinds.append('N')
    for m in V.monomials(degQ, pA):
        cols.append(V.shift(Fq1, m, 1, p)); kinds.append('Q')
    nNQ = len(cols)
    D2 = V.mul(Dp, Dp, p)
    if with_target_cols:
        for m in V.monomials(6*q, 1):
            cols.append(V.shift(D2, m, 1, p)); kinds.append('A')
    rhs = V.scale(V.mul(Tp, D2, p), p - 1, p) if Tp else {}
    rows = {}
    for c in cols:
        for e in c: rows.setdefault(e, len(rows))
    for e in rhs: rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    return cols, kinds, rows, rhs, nNQ, R, C

def dense(cols, rows, rhs, R, C):
    flat = [0]*(R*(C + 1))
    W = C + 1
    for j, c in enumerate(cols):
        for e, v in c.items(): flat[rows[e]*W + j] = v
    for e, v in rhs.items(): flat[rows[e]*W + C] = v
    return nmod_mat(R, C + 1, flat, p)

def consistent(Fx, D, target, label):
    t0 = time.time()
    cols, kinds, rows, rhs, nNQ, R, C = build(Fx, D, target)
    Maug = dense(cols, rows, rhs, R, C)
    r_aug = Maug.rank()
    M0 = dense(cols, rows, {}, R, C)   # last column zero
    r0 = M0.rank()
    out = dict(label=label, d=V.deg(D), rows=R, cols=C, rank=r0, rank_aug=r_aug, consistent=(r0 == r_aug), nullity=C - r0,
               p=p, seconds=round(time.time() - t0, 1))
    with open('/tmp/hodge/results/quadratic-cubic-pullback/verify/fixedD_own.jsonl', 'a') as f: f.write(json.dumps(out) + '\n')
    return out

def served_dim(Fx, D, label):
    t0 = time.time()
    cols, kinds, rows, rhs, nNQ, R, C = build(Fx, D, None, with_target_cols=True)
    M = dense(cols, rows, {}, R, C)
    rM = M.rank()
    M0 = dense(cols[:nNQ], rows, {}, R, nNQ)
    r0 = M0.rank()
    nA = C - nNQ
    sd = nA - (rM - r0)
    # is the specific target A in the served space? (consistency of the inhomogeneous system)
    out = dict(label=label, d=V.deg(D), rows=R, cols_NQ=nNQ, cols_A=nA, rank_NQ=r0, rank_all=rM, served_dim=sd, p=p,
               seconds=round(time.time() - t0, 1))
    with open('/tmp/hodge/results/quadratic-cubic-pullback/verify/fixedD_own.jsonl', 'a') as f: f.write(json.dumps(out) + '\n')
    return out

def P(s): return V.parse(s)
Fi_ex = [V.pderiv(F, i) for i in range(6)]
def lin(cs, Fs):
    out = {}
    for c, Fx in zip(cs, Fs): out = V.add(out, V.mul(P(c) if isinstance(c, str) else c, Fx))
    return out
join012 = V.add(V.add(V.power(Q[0], 3), V.power(Q[1], 3)), V.power(Q[2], 3))
join345 = V.add(V.add(V.power(Q[3], 3), V.power(Q[4], 3)), V.power(Q[5], 3))
case = sys.argv[1]
if case == 'fermat':
    o = consistent(F0, P('X0**6+X1**6+X2**6'), A, 'fermat control D=X0^6+X1^6+X2^6 target A')
    log('C14 Fermat control: fixed-D system (own builder) consistent for D=X0^6+X1^6+X2^6, target A', 'CONFIRMED' if o['consistent'] else 'REFUTED', json.dumps(o))
    o = served_dim(F0, P('X0**6+X1**6+X2**6'), 'fermat control served dim')
    log('C15 Fermat control: served dimension 118 = 18 + 100 for D=X0^6+X1^6+X2^6 (own builder)', 'CONFIRMED' if o['served_dim'] == 118 else 'REFUTED', json.dumps(o))
elif case == 'join012':
    o = consistent(F, join012, A, 'join012 target A')
    log('C16 join012 = Q0^3+Q1^3+Q2^3 does not serve A (fixed-D system inconsistent; own builder, second prime)', 'CONFIRMED' if not o['consistent'] else 'REFUTED', json.dumps(o))
    o = consistent(F, join012, delta, 'join012 target delta')
    log('C17 join012 serves delta: fixed-D system with target delta consistent (own builder, second prime; mod-p evidence, not an exact lift)', 'CONFIRMED' if o['consistent'] else 'REFUTED', json.dumps(o))
    o = served_dim(F, join012, 'join012 served dim')
    log('C18 join012 served dimension = 19 = 18 + 1 (own builder, second prime)', 'CONFIRMED' if o['served_dim'] == 19 else 'REFUTED', json.dumps(o))
elif case == 'join345':
    o = consistent(F, join345, delta, 'join345 target delta')
    log('C19 join345 serves delta (own builder, second prime; mod-p evidence)', 'CONFIRMED' if o['consistent'] else 'REFUTED', json.dumps(o))
    o = served_dim(F, join345, 'join345 served dim')
    log('C19b join345 served dimension (own builder)', 'CONFIRMED' if o['served_dim'] == 19 else 'WEAKENED', json.dumps(o))
elif case == 'polar5':
    for lab, D in [('F_0', Fi_ex[0]), ('F_0+2F_1+3F_2', lin(['1', '2', '3'], Fi_ex[:3])), ('F_3 (Y-odd)', Fi_ex[3]),
                   ('F_1-F_2', lin(['1', '-1'], [Fi_ex[1], Fi_ex[2]])), ('F_4+F_5 (Y-odd)', lin(['1', '1'], [Fi_ex[4], Fi_ex[5]]))]:
        o = consistent(F, D, A, f'polar5 {lab} target A')
        o2 = served_dim(F, D, f'polar5 {lab} served dim')
        log(f'C20 degree-5 polar D={lab}: A not served, served dim 18 (own builder, second prime; mod-p evidence)',
            'CONFIRMED' if (not o['consistent'] and o2['served_dim'] == 18) else 'REFUTED', json.dumps(o) + ' ' + json.dumps(o2))
elif case == 'polar6':
    import random
    random.seed(20261102)
    def rnd_lin(par):
        idx = range(0, 3) if par == 0 else range(3, 6)
        return {tuple(1 if k == i else 0 for k in range(6)): V.Fraction(random.randint(-9, 9)) for i in idx}
    Dr = {}
    for i in range(6):
        Dr = V.add(Dr, V.mul(rnd_lin(0 if i < 3 else 1), Fi_ex[i]))
    for lab, D in [('X0F0 - Y0F3', V.add(V.mul(P('X0'), Fi_ex[0]), V.scale(V.mul(P('Y0'), Fi_ex[3]), -1))), ('random Y-even J6 element', Dr),
                   ('X1F0 + X0F1 + Y1F3 + Y0F4', lin(['X1', 'X0', 'Y1', 'Y0'], [Fi_ex[0], Fi_ex[1], Fi_ex[3], Fi_ex[4]]))]:
        o = consistent(F, D, A, f'polar6 {lab} target A')
        o2 = served_dim(F, D, f'polar6 {lab} served dim')
        log(f'C21 degree-6 polar D={lab}: A not served, served dim 18 (own builder, second prime; mod-p evidence)',
            'CONFIRMED' if (not o['consistent'] and o2['served_dim'] == 18) else 'REFUTED', json.dumps(o) + ' ' + json.dumps(o2))
elif case == 'delta':
    o = consistent(F, delta, A, 'delta target A')
    log('C22 D=delta does not serve A (own builder, second prime)', 'CONFIRMED' if not o['consistent'] else 'REFUTED', json.dumps(o))
elif case == 'polar7':
    D = lin(['Y0**2', 'X0*Y0'], [Fi_ex[0], Fi_ex[3]])
    o = consistent(F, D, A, 'polar7 Y0^2F0 + X0Y0F3 target A')
    log('C23 degree-7 polar D=Y0^2F0+X0Y0F3: A not served (own builder, second prime; mod-p evidence)', 'CONFIRMED' if not o['consistent'] else 'REFUTED', json.dumps(o))
print('done', case)
