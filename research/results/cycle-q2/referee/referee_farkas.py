#!/usr/bin/env python3
"""referee_farkas.py -- INDEPENDENT exact replay of stored Farkas certificates (lane cycle-q2, wave2-excl/certs/*.json).

No engine, no sector restriction, no modular arithmetic.  For a certificate (D, row monomials r, rational y_r) we check with
exact integer arithmetic that the linear functional Y(P) = sum_r y_r [P]_r annihilates EVERY column of the full fixed-D system
    L(N, Q) = -q D sum_i N_i F2_i + F2 (D sum_i d_i N_i - sum_i N_i d_i D) - F2^(q+1) Q       (F2 = 2F, integer coefficients)
over all monomials N_i of degree deg D + 6q - 5 and all monomials Q of degree 2 deg D - 6, and that Y(A D^2) != 0.
Since the certificate identity is  A D^2 + L(N, Q) = 0, this proves that no (N, Q) over any field of characteristic 0 exists
for this D (scaling F by 2 only rescales columns).  The column functionals are computed by correlation (adjoint of
multiplication), so the full column set is covered without enumerating it.
usage: python3 -I referee_farkas.py certs/B6.json [more.json ...]   (writes referee_farkas.results.jsonl)
"""
import sys, json, ast, time, math, hashlib
from fractions import Fraction
sys.set_int_max_str_digits(0)
VARS = ['X0','X1','X2','Y0','Y1','Y2']
Z6 = (0,)*6

def padd(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + s*v
    return {k: v for k, v in r.items() if v}
def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(x+y for x, y in zip(k1, k2)); r[k] = r.get(k, 0) + v1*v2
    return {k: v for k, v in r.items() if v}
def pscale(a, c): return {k: v*c for k, v in a.items() if v*c}
def pderiv(P, i):
    r = {}
    for k, v in P.items():
        if k[i]:
            e = list(k); e[i] -= 1; r[tuple(e)] = r.get(tuple(e), 0) + v*k[i]
    return {k: v for k, v in r.items() if v}
def ppow(P, n):
    r = {Z6: Fraction(1)}
    for _ in range(n): r = pmul(r, P)
    return r
def parse(s):
    def ev(n):
        if isinstance(n, ast.Expression): return ev(n.body)
        if isinstance(n, ast.Constant): return {Z6: Fraction(n.value)}
        if isinstance(n, ast.Name): e = [0]*6; e[VARS.index(n.id)] = 1; return {tuple(e): Fraction(1)}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return pscale(ev(n.operand), -1)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.UAdd): return ev(n.operand)
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                assert isinstance(n.right, ast.Constant); return ppow(ev(n.left), int(n.right.value))
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return padd(a, b)
            if isinstance(n.op, ast.Sub): return padd(a, b, -1)
            if isinstance(n.op, ast.Mult): return pmul(a, b)
            if isinstance(n.op, ast.Div):
                assert list(b) == [Z6]; return pscale(a, 1/b[Z6])
        raise ValueError(ast.dump(n))
    return ev(ast.parse(s.replace('^', '**'), mode='eval'))
def degree(P): return max(sum(k) for k in P)

def corr(P, Y):
    """adjoint of multiplication: corr(P,Y)[m] = sum_s P[s] * Y[s+m]  (= Y(P * m))."""
    out = {}
    for r, yv in Y.items():
        for s, pv in P.items():
            if all(a >= b for a, b in zip(r, s)):
                m = tuple(a-b for a, b in zip(r, s)); out[m] = out.get(m, 0) + pv*yv
    return out

def check(path, F, A, q):
    t0 = time.time()
    c = json.load(open(path))
    assert c['F'] == '2F' and c['q'] == q
    D = parse(c['D']); d = degree(D)
    assert all(v.denominator == 1 for v in D.values()), 'D not integral'
    D = {k: int(v) for k, v in D.items()}
    if c.get('y') is None:
        return {'cert': path, 'status': 'NO-VECTOR-STORED (y omitted)', 'd': d}
    rows = [tuple(r) for r in c['row_monomials']]
    y = [Fraction(v) for v in c['y']]
    assert len(rows) == len(y)
    L = 1
    for v in y: L = L*v.denominator//math.gcd(L, v.denominator)
    Y = {r: int(v*L) for r, v in zip(rows, y) if v != 0}          # integral scaling of y (does not affect the certificate)
    assert all(sum(r) == 2*d + 6*q for r in Y), 'row degree'
    F2 = {k: int(2*v) for k, v in F.items()}
    assert all(v.denominator == 1 for v in (Fraction(2*v) for v in F.values()))
    A_ = {k: int(v) for k, v in A.items()}
    F2i = [pderiv(F2, i) for i in range(6)]
    Di = [pderiv(D, i) for i in range(6)]
    F2D = pmul(F2, D)
    T2 = corr(F2D, Y)                       # indexed by m' (deg d+6q-6):  Y(F2 D m')
    bad = 0; ncols = 0
    for i in range(6):
        T1 = corr(pmul(D, F2i[i]), Y)       # Y(D F2_i m)
        T3 = corr(pmul(F2, Di[i]), Y)       # Y(F2 d_iD m)
        keys = set(T1) | set(T3)
        for mp in T2:
            e = list(mp); e[i] += 1; keys.add(tuple(e))
        for m in keys:
            val = -q*T1.get(m, 0) - T3.get(m, 0)
            if m[i]:
                e = list(m); e[i] -= 1; val += m[i]*T2.get(tuple(e), 0)
            ncols += 1
            if val != 0: bad += 1
    TQ = corr(ppow({k: Fraction(v) for k, v in F2.items()}, q+1), Y)
    badQ = sum(1 for v in TQ.values() if v != 0); ncols += len(TQ)
    AD2 = pmul(A_, pmul(D, D))
    yb = sum(Y.get(r, 0)*v for r, v in AD2.items())
    ok = (bad == 0 and badQ == 0 and yb != 0)
    h = hashlib.sha256(('\n'.join(c['y']) + '\n' + json.dumps(c['row_monomials'])).encode()).hexdigest()
    return {'cert': path, 'status': 'VERIFIED' if ok else 'FAILED', 'd': d, 'nrows_y': len(Y), 'ncols_touched': ncols,
            'bad_N_cols': bad, 'bad_Q_cols': badQ, 'y_dot_AD2_nonzero': yb != 0, 'y_dot_AD2_bits': yb.bit_length(),
            'sha256_matches_file': h == c['sha256_y'], 'D': c['D'][:120], 'seconds': round(time.time()-t0, 1)}

if __name__ == '__main__':
    rec = open('/home/kejace/devel/sextics/data/cycle-q2/instance.rec').read()
    F = parse(rec.split('\nF: ')[1].split('\n')[0]); A = parse(rec.split('\nA: ')[1].split('\n')[0]); q = int(rec.split('\nq: ')[1].split('\n')[0])
    for p in sys.argv[1:]:
        r = check(p, F, A, q)
        print(json.dumps(r), flush=True)
        with open('/tmp/hodge/referee/cycle-q2/referee_farkas.results.jsonl', 'a') as f: f.write(json.dumps(r) + '\n')
