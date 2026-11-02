"""Exact (mod p) normalized E2 test for single denominators (q = 1).

  python3 wave2/e2exact.py --expr "<D expression>" [--tag e2exact] [--name label]
  python3 wave2/e2exact.py --kernel-json out/e1kernel_d7_p0.json --index 0      # kernel representative mod p
  python3 wave2/e2exact.py --members-from out/structured.jsonl --family "vars[X0X1X2Y0]" --d 10 --parity 0

Steps: E1 lift N0 (augmented rref of the Jacobian generators in degree d+6), then the linear system
   T_D(c) + F C = -B_D(N0)   over F_p   (columns: Koszul c_ij of degree d-4 with parity, C in S_{2d-6}^{pA}),
consistency / rank / nullity by nmod_mat rref (dense; optional random row compression when rows >> cols).
"""
import sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
from flint import nmod_mat
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--expr', default=None)
ap.add_argument('--kernel-json', default=None)
ap.add_argument('--index', type=int, default=0)
ap.add_argument('--members-from', default=None)
ap.add_argument('--family', default=None)
ap.add_argument('--d', type=int, default=None)
ap.add_argument('--parity', type=int, default=None)
ap.add_argument('--prime-index', type=int, default=0)
ap.add_argument('--tag', default='e2exact')
ap.add_argument('--name', default='')
ap.add_argument('--compress', action='store_true', help='random row compression to ncols+40 rows before rref')
ap.add_argument('--random-member', action='store_true', help='D = random mod-p combination of J-type generators and kernel reps of degree --d, parity --parity')
ap.add_argument('--seed', type=int, default=0)
args = ap.parse_args()

F, A, q = load()
p = core.PRIMES[args.prime_index]
Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
pA = core.parity_of(A)


def e1_lift(D):
    d = core.degree(D); pD = core.parity_of(D)
    pT = (pA + pD) % 2
    gens = jacobian.jacobian_generators(F, d + 6, p, pT)
    rhs = mul(Ap, D, p)
    rows = {}
    for c in [g[2] for g in gens] + [rhs]:
        for e in c:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(gens)
    flat = [0] * (R * (C + 1))
    for j, g in enumerate(gens):
        for e, v in g[2].items():
            flat[rows[e] * (C + 1) + j] = v
    for e, v in rhs.items():
        flat[rows[e] * (C + 1) + C] = v
    Mr, rank = nmod_mat(R, C + 1, flat, p).rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    if pivots and pivots[-1] == C:
        return None, {'e1': False, 'rank': rank, 'rows': R, 'cols': C}
    N = [{} for _ in range(NV)]
    for i, pv in enumerate(pivots):
        v = int(Mr[i, C])
        if v:
            gi, gm, _ = gens[pv]
            N[gi][gm] = (N[gi].get(gm, 0) + v) % p
    NF = {}
    for i in range(NV):
        NF = add(NF, mul(N[i], Fi[i], p), p)
    assert sub(NF, rhs, p) == {}, 'lift verification failed'
    return N, {'e1': True, 'rank': rank, 'rows': R, 'cols': C}


def e2_system(D, N0):
    d = core.degree(D); pD = core.parity_of(D)
    Srow = Space(2 * d, pA)
    SC = Space(2 * d - 6, pA)
    Di = [pderiv(D, i, p) for i in range(NV)]
    DF = [mul(D, Fi[j], p) for j in range(NV)]
    cols = []
    labels = []
    for i in range(NV):
        for j in range(i + 1, NV):
            pc = (pA - pD + par_idx(i) + par_idx(j)) % 2
            Pij = sub(mul(Fi[j], Di[i], p), mul(Fi[i], Di[j], p), p)
            for m in monomials(d - 4, pc):
                col = {}
                if m[i]:
                    m2 = list(m); m2[i] -= 1
                    col = add(col, shift(DF[j], tuple(m2), m[i], p), p)
                if m[j]:
                    m2 = list(m); m2[j] -= 1
                    col = sub(col, shift(DF[i], tuple(m2), m[j], p), p)
                col = sub(col, shift(Pij, m, 1, p), p)
                cols.append(col); labels.append(('c', i, j, m))
    nc = len(cols)
    for m in monomials(2 * d - 6, pA):
        cols.append(shift(Fp, m, 1, p)); labels.append(('C', m))
    rhs = scale(B_of(D, N0, p), p - 1, p)
    return Srow, cols, labels, nc, rhs


def solve(D, name):
    t0 = time.time()
    d = core.degree(D); pD = core.parity_of(D)
    N0, e1info = e1_lift(D)
    rec = {'batch': 'e2exact', 'name': name, 'd': d, 'pD': pD, 'p': p, 'D': core.to_str(jacobian.symmetric(D, p))[:2000], 'nterms_D': len(D), 'e1': e1info}
    if N0 is None:
        rec.update({'consistent': False, 'note': 'e1 fails', 'time_s': round(time.time() - t0, 1)})
        return rec
    t1 = time.time()
    Srow, cols, labels, nc, rhs = e2_system(D, N0)
    rows = {}
    for c in cols + [rhs]:
        for e in c:
            rows.setdefault(e, len(rows))
    R, C = len(rows), len(cols)
    rec.update({'rows': R, 'cols': C, 'n_c': nc, 'n_C': C - nc, 't_build': round(time.time() - t1, 1)})
    t2 = time.time()
    if args.compress and R > C + 60:
        import scipy.sparse as sp
        ri, ci, vals = [], [], []
        for j, c in enumerate(cols):
            for e, v in c.items():
                ri.append(rows[e]); ci.append(j); vals.append(v)
        for e, v in rhs.items():
            ri.append(rows[e]); ci.append(C); vals.append(v)
        Msp = sp.csr_matrix((np.array(vals, dtype=np.int64), (np.array(ri), np.array(ci))), shape=(R, C + 1))
        k = C + 40
        rng = np.random.default_rng(1)
        Rm = rng.integers(0, p, size=(k, R), dtype=np.int64)
        MT = Msp.T.tocsr()
        lo = (MT @ (Rm & 0xFFFF).T) % p
        hi = (MT @ (Rm >> 16).T) % p
        Cm = ((lo + (hi * (1 << 16)) % p) % p).T
        Mn = nmod_mat(k, C + 1, np.ascontiguousarray(Cm).ravel().tolist(), p)
        rec['compressed_rows'] = k
    else:
        flat = [0] * (R * (C + 1))
        for j, c in enumerate(cols):
            for e, v in c.items():
                flat[rows[e] * (C + 1) + j] = v
        for e, v in rhs.items():
            flat[rows[e] * (C + 1) + C] = v
        Mn = nmod_mat(R, C + 1, flat, p)
    Mr, rank = Mn.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= C and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    # rank of the coefficient matrix alone
    rank0 = rank - (0 if consistent else 1)
    rec.update({'consistent': bool(consistent), 'rank': int(rank0), 'rank_aug': int(rank), 'nullity': int(C - rank0),
                't_rref': round(time.time() - t2, 1), 'time_s': round(time.time() - t0, 1)})
    return rec


if __name__ == '__main__':
    Ds = []
    if args.expr:
        Ds.append((args.name or args.expr[:40], core.reduce_mod(core.parse(args.expr), p)))
    if args.kernel_json:
        kj = json.load(open(args.kernel_json)); assert kj['p'] == p
        b = kj['basis_modp'][args.index]
        Ds.append((args.name or f'{os.path.basename(args.kernel_json)}#{args.index}', core.parse(b, p)))
    if args.members_from:
        for line in open(args.members_from):
            o = json.loads(line)
            if o.get('family') == args.family and o.get('d') == args.d and o.get('parity') == args.parity:
                for k, m in enumerate(o.get('members_modp', [])):
                    Ds.append((f'{args.family} d={args.d} par={args.parity} #{k}', core.parse(m, p)))
    if args.random_member:
        import random
        rnd = random.Random(args.seed)
        d, par = args.d, args.parity
        D = {}
        for i in range(NV):
            for m in monomials(d - 5, (par - par_idx(i)) % 2):
                D = add(D, shift(Fi[i], m, rnd.randrange(p), p), p)
        kpath = f'{OUT}/e1kernel_d{d}_p{par}.json'
        nker = 0
        if os.path.exists(kpath):
            kj = json.load(open(kpath)); assert kj['p'] == p
            for b in kj['basis_modp']:
                D = add(D, scale(core.parse(b, p), rnd.randrange(p), p), p); nker += 1
        Ds.append((f'random member d={d} par={par} seed={args.seed} (J + {nker} kernel)', D))
    for name, D in Ds:
        rec = solve(D, name)
        print(json.dumps(rec), flush=True)
        log_line(f'{OUT}/{args.tag}.jsonl', rec)
