"""Exp C: formal continuation of a Fermat certificate along the pencil F_c = Fermat + c*P (P = 7/5 * coupling),
A = X0^3 Y0^3, q = 1, deg D = 6.  Order-k equation: L(D_k,N_k,Q_k) = -(known lower-order terms), where L is the
linearization of E = A D^2 - D N(F) + F B_D(N) - F^2 Q at the base certificate (D_0,N_0,Q_0) w.r.t. F_0.
Klein-sector restricted (D_k: Y-even,Y2-even; N_{k,i}: parity of A*Z_i; Q_k: parity of A), modulo a prime p.
At order k+1 the free (kernel) part of order k is also allowed to vary (linear).
Usage: python3 fermat_cont.py <name> <prime> <base> <maxorder>     base in {X, Y}
"""
import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint
name, p, base, maxorder = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
t0 = time.time()
F6, A, q, rec = core.load_instance('reflection-coupled-residual')
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
Pq = core.sub(F6, F0)                      # c*P with c = 5/7  -> we use the pencil F_0 + t*(F6 - F0); t=1 is the instance
Pm = core.reduce_mod(Pq, p); F0m = core.reduce_mod(F0, p); Am = core.reduce_mod(A, p)
def mul(a, b): return core.mul(a, b, p)
def add(a, b): return core.add(a, b, p)
def sub(a, b): return core.sub(a, b, p)
def scale(a, c): return core.scale(a, c % p, p)
def pd(a, i): return core.pderiv(a, i, p)
def apply(N, G):  # N(G) = sum N_i d_i G
    out = {}
    for i in range(6):
        if N[i]: out = add(out, mul(N[i], pd(G, i)))
    return out
def div(N):
    out = {}
    for i in range(6): out = add(out, pd(N[i], i))
    return out
def B(D, N):  # D div N - N(D)
    return sub(mul(D, div(N)), apply(N, D))
inv6 = pow(6, -1, p)
if base == 'X':
    I = [0, 1, 2]
else:
    I = [3, 4, 5]
D0 = {}
for i in I:
    e = [0]*6; e[i] = 6; D0[tuple(e)] = 1
N0 = [{} for _ in range(6)]
for i in I:
    e = [0]*6; e[i] = 1
    N0[i] = scale(core.mul(Am, {tuple(e): 1}), inv6)
Q0 = {}
# check base identity
E0 = sub(add(sub(mul(Am, mul(D0, D0)), mul(D0, apply(N0, F0m))), mul(F0m, B(D0, N0))), mul(mul(F0m, F0m), Q0))
assert not E0, "base certificate fails"
print(f'[{name}] base {base} ok', flush=True)
# parity helpers
def yp(e): return (e[3]+e[4]+e[5]) % 2
def y2p(e): return e[5] % 2
parA = (1, 0)   # (Y-parity, Y2-parity) of A
def par_of(poly):
    e = next(iter(poly)); return (yp(e), y2p(e))
# unknown monomial bases
D_mons = [e for e in core.monomials(6) if (yp(e), y2p(e)) == (0, 0)]
N_mons = []
for i in range(6):
    ei = [0]*6; ei[i] = 1
    want = ((parA[0] + yp(ei)) % 2, (parA[1] + y2p(ei)) % 2)
    N_mons.append([e for e in core.monomials(7) if (yp(e), y2p(e)) == want])
Q_mons = [e for e in core.monomials(6) if (yp(e), y2p(e)) == parA]
row_mons = [e for e in core.monomials(18) if (yp(e), y2p(e)) == parA]
row_index = {e: k for k, e in enumerate(row_mons)}
cols = [('D', None, e) for e in D_mons] + [('N', i, e) for i in range(6) for e in N_mons[i]] + [('Q', None, e) for e in Q_mons]
print(f'[{name}] rows={len(row_mons)} cols={len(cols)} (D {len(D_mons)}, N {sum(len(m) for m in N_mons)}, Q {len(Q_mons)})', flush=True)
N0F0 = apply(N0, F0m)           # = A D0
F0sq = mul(F0m, F0m)
def L_apply(D1, N1, Q1):
    # 2 A D0 D1 - D1 N0(F0) - D0 N1(F0) + F0 (B_{D0}(N1) + B_{D1}(N0)) - F0^2 Q1
    t = scale(mul(Am, mul(D0, D1)), 2)
    t = sub(t, mul(D1, N0F0))
    t = sub(t, mul(D0, apply(N1, F0m)))
    t = add(t, mul(F0m, add(B(D0, N1), B(D1, N0))))
    t = sub(t, mul(F0sq, Q1))
    return t
# build matrix of L (columns)
tL = time.time()
M = flint.nmod_mat(len(row_mons), len(cols), p)
for j, (kind, i, e) in enumerate(cols):
    if kind == 'D':
        col = L_apply({e: 1}, [{}]*6, {})
    elif kind == 'N':
        N1 = [{}]*6; N1[i] = {e: 1}
        col = L_apply({}, N1, {})
    else:
        col = L_apply({}, [{}]*6, {e: 1})
    for m, cf in col.items():
        M[row_index[m], j] = cf
print(f'[{name}] L built in {time.time()-tL:.1f}s', flush=True)
rankL = M.rank()
print(f'[{name}] rank L = {rankL}, nullity = {len(cols)-rankL}', flush=True)
def vec(poly):
    v = flint.nmod_mat(len(row_mons), 1, p)
    for m, cf in poly.items():
        if m not in row_index:
            raise ValueError('parity violation in RHS: %s' % (m,))
        v[row_index[m], 0] = cf
    return v
def unpack(x):
    D = {}; N = [{} for _ in range(6)]; Q = {}
    for j, (kind, i, e) in enumerate(cols):
        c = int(x[j, 0]) if hasattr(x, 'nrows') else int(x[j])
        if not c: continue
        if kind == 'D': D[e] = c
        elif kind == 'N': N[i][e] = c
        else: Q[e] = c
    return D, N, Q
def solve(Mat, rhs):
    """Return (consistent, particular solution as list) for Mat x = rhs mod p via rref of augmented matrix."""
    nr, nc = Mat.nrows(), Mat.ncols()
    Aug = flint.nmod_mat(nr, nc + 1, p)
    for r in range(nr):
        for c in range(nc):
            Aug[r, c] = Mat[r, c]
        Aug[r, nc] = rhs[r, 0]
    R, rk = Aug.rref()
    # inconsistent iff a pivot in the last column
    pivots = []
    for r in range(rk):
        c = 0
        while c <= nc and R[r, c] == 0: c += 1
        pivots.append(c)
    if pivots and pivots[-1] == nc:
        return False, None, rk
    x = [0]*nc
    for r, c in enumerate(pivots):
        x[c] = int(R[r, nc])
    return True, x, rk
# orders
Ds = [D0]; Ns = [N0]; Qs = [Q0]
results = {'name': name, 'p': p, 'base': base, 'rankL': rankL, 'nullity': len(cols) - rankL, 'orders': []}
for k in range(1, maxorder + 1):
    tk = time.time()
    # known terms at order k excluding the L part
    rest = {}
    for i in range(1, k):      # A D_i D_j with i+j=k, 0<i,j<k
        j = k - i
        rest = add(rest, mul(Am, mul(Ds[i], Ds[j])))
    for i in range(1, k):      # - D_i N_j(F0), i+j=k, 0<i,j<k
        j = k - i
        rest = sub(rest, mul(Ds[i], apply(Ns[j], F0m)))
    for i in range(0, k):      # - D_i N_j(P), i+j=k-1
        j = k - 1 - i
        rest = sub(rest, mul(Ds[i], apply(Ns[j], Pm)))
    for i in range(1, k):      # F0 B_{D_i}(N_j), i+j=k, 0<i,j<k
        j = k - i
        rest = add(rest, mul(F0m, B(Ds[i], Ns[j])))
    for i in range(0, k):      # P B_{D_i}(N_j), i+j=k-1
        j = k - 1 - i
        rest = add(rest, mul(Pm, B(Ds[i], Ns[j])))
    if k - 1 >= 0: rest = sub(rest, scale(mul(mul(F0m, Pm), Qs[k-1]), 2))
    if k - 2 >= 0: rest = sub(rest, mul(mul(Pm, Pm), Qs[k-2]))
    rhs = vec(scale(rest, p - 1))     # L x = -rest
    ok, x, rk = solve(M, rhs)
    info = {'order': k, 'consistent': ok, 'rank_aug': rk, 'seconds': round(time.time() - tk, 1)}
    print(f'[{name}] order {k}: consistent={ok} rank_aug={rk} ({time.time()-tk:.1f}s)', flush=True)
    results['orders'].append(info)
    json.dump(results, open(f'{name}.json', 'w'))
    if not ok:
        break
    Dk, Nk, Qk = unpack(x)
    Ds.append(Dk); Ns.append(Nk); Qs.append(Qk)
    results['orders'][-1]['nnz'] = [len(Dk), sum(len(n) for n in Nk), len(Qk)]
results['seconds'] = time.time() - t0
json.dump(results, open(f'{name}.json', 'w'))
print(json.dumps(results), flush=True)
