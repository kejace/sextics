"""Bundle certificate: A = sum_S A_S with A_S served by D_S; exact (N_S,Q_S) via joint modular solve + CRT.

usage: python3 bundle.py <served_rows_tag.pkl> [maxk]
Uses the pickled served rows (mod p) to pick a small subset of denominators whose served spaces
span A, then solves the joint linear system exactly and writes witness-candidate.json (combined)
and bundle-pieces.json (each piece).
"""
import sys, time, json, pickle, os, random
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
from fractions import Fraction
from flint import nmod_mat
import numpy as np
import scipy.sparse as sp

WD = qcp.WD
F_, A_, q_ = qcp.F, qcp.A, 1
delta = P(open(f'{WD}/delta.txt').read())


def Fpart(xdeg):
    return {e: c for e, c in F_.items() if e[0] + e[1] + e[2] == xdeg}


def EI(I, G=None):
    G = F_ if G is None else G
    out = {}
    for i in I:
        out = core.add(out, core.mul(qcp.var(i), core.pderiv(G, i)))
    return out


env = {'P': P, 'Q': qcp.Q, 'core': core, 'delta': delta, 'F': F_, 'Fpart': Fpart, 'EI': EI}

pk = pickle.load(open(sys.argv[1], 'rb'))
maxk = int(sys.argv[2]) if len(sys.argv) > 2 else 6
p0 = pk['p']; odd6 = pk['odd6']; rows = pk['rows']; owners = pk['owners']; items = pk['items']
col = {e: j for j, e in enumerate(odd6)}
Avec = [0] * len(odd6)
for e, v in core.reduce_mod(A_, p0).items():
    Avec[col[e]] = v
names = []
Dexprs = {}
for line in items:
    parts = line.split('|')
    names.append(parts[0].strip()); Dexprs[parts[0].strip()] = parts[1].strip()


def span_contains(rowsubset, vec, p):
    if not rowsubset:
        return False
    M = nmod_mat(len(rowsubset), len(odd6), [x for r in rowsubset for x in r], p)
    M2 = nmod_mat(len(rowsubset) + 1, len(odd6), [x for r in rowsubset for x in r] + vec, p)
    return M.rank() == M2.rank()


# greedy selection: add denominators (in given order) until A is in the span; then prune
chosen = []
for nm in names:
    chosen.append(nm)
    sub = [r for r, o in zip(rows, owners) if o in chosen]
    if span_contains(sub, Avec, p0):
        break
else:
    print('A not in the union of served spaces; abort'); sys.exit(1)
# prune
for nm in list(chosen):
    trial = [c for c in chosen if c != nm]
    sub = [r for r, o in zip(rows, owners) if o in trial]
    if span_contains(sub, Avec, p0):
        chosen = trial
print('chosen denominators:', chosen, flush=True)
if len(chosen) > maxk:
    print('too many pieces'); sys.exit(1)
Ds = [eval(Dexprs[nm], env) for nm in chosen]
Ds = [P(D) if isinstance(D, str) else D for D in Ds]
k = len(Ds)


def joint_solve(p, seed=1):
    """Build joint compressed system mod p; return dict with consistent, pivots, x, layout."""
    blocks = []
    for D in Ds:
        S = fixed_d.FixedDSystem(F_, A_, q_, D, p, with_target_space=True)
        blocks.append(S)
    # per block: compress rows to (ncols + 40)
    rng = np.random.default_rng(seed)
    mats = []
    layout = []  # (block index, S.cols) ; global column offsets
    off = 0
    for bi, S in enumerate(blocks):
        M = fixed_d._sparse_matrix(S, with_rhs=False)   # rows x cols (cols = N,Q,A columns)
        kk = min(S.nrows, S.ncols + 40)
        Cm = fixed_d._compress(M, p, kk, seed + bi)       # kk x ncols
        mats.append(Cm)
        layout.append((off, S.ncols))
        off += S.ncols
    totcols = off
    # coupling rows: sum_S A_S = A  (226 rows), A columns in each block are those with kind 'A'
    nA = len(odd6)
    R = sum(m.shape[0] for m in mats) + nA
    flat = np.zeros((R, totcols + 1), dtype=np.int64)
    r0 = 0
    for bi, (Cm, (off_, nc)) in enumerate(zip(mats, layout)):
        flat[r0:r0 + Cm.shape[0], off_:off_ + nc] = Cm
        r0 += Cm.shape[0]
    aidx = {e: i for i, e in enumerate(odd6)}
    for bi, S in enumerate(blocks):
        off_, nc = layout[bi]
        for j, c in enumerate(S.cols):
            if c[0] == 'A':
                flat[r0 + aidx[c[2]], off_ + j] = 1
    for e, v in core.reduce_mod(A_, p).items():
        flat[r0 + aidx[e], totcols] = v
    M = nmod_mat(R, totcols + 1, flat.ravel().tolist(), p)
    Mr, rank = M.rref()
    pivots = []
    for i in range(rank):
        j = 0
        while j <= totcols and int(Mr[i, j]) == 0:
            j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == totcols)
    x = None
    if consistent:
        x = [0] * totcols
        for i, j in enumerate(pivots):
            x[j] = int(Mr[i, totcols])
    return {'consistent': consistent, 'rank': rank, 'pivots': pivots, 'x': x, 'blocks': blocks,
            'layout': layout, 'totcols': totcols}


def unpack(blocks, layout, vals):
    pieces = []
    for bi, S in enumerate(blocks):
        off_, nc = layout[bi]
        N = [{} for _ in range(6)]; Qp = {}; AS = {}
        for j, c in enumerate(S.cols):
            v = vals[off_ + j]
            if v == 0:
                continue
            if c[0] == 'N': N[c[1]][c[2]] = v
            elif c[0] == 'Q': Qp[c[2]] = v
            else: AS[c[2]] = v
        pieces.append((N, Qp, AS))
    return pieces


t0 = time.time()
acc = None; piv0 = None
for kp, p in enumerate(core.PRIMES[:10]):
    r = joint_solve(p)
    print(f'[p={p}] consistent={r["consistent"]} rank={r["rank"]} cols={r["totcols"]} ({time.time()-t0:.0f}s)', flush=True)
    if not r['consistent']:
        print('joint system inconsistent mod p -- served-space analysis disagrees?'); sys.exit(1)
    if piv0 is None or piv0 != r['pivots']:
        if piv0 is not None:
            print('pivot structure changed; restarting CRT', flush=True)
        piv0 = r['pivots']; acc = (r['x'], p)
    else:
        xs, m = acc
        xs = [core.crt_pair(a, m, b, p)[0] for a, b in zip(xs, r['x'])]
        acc = (xs, m * p)
    xs, m = acc
    rat = [core.ratrecon(a, m) for a in xs]
    if all(v is not None for v in rat):
        pieces = unpack(r['blocks'], r['layout'], rat)
        ok_all = True
        for (N, Qp, AS), D in zip(pieces, Ds):
            ok, info = core.verify_witness(F_, AS, 1, D, N, Qp)
            print(f'   piece d={core.degree(D)}: verify={ok} {info}', flush=True)
            ok_all &= ok
        Asum = {}
        for (_, _, AS) in pieces:
            Asum = core.add(Asum, AS)
        ok_sum = core.sub(Asum, A_) == {}
        print(f'   sum of pieces == A: {ok_sum}', flush=True)
        if ok_all and ok_sum:
            # combine via addition formula
            Dt, Nt, Qt = Ds[0], pieces[0][0], pieces[0][1]
            for (N, Qp, AS), D in list(zip(pieces, Ds))[1:]:
                # (Dt,Nt,Qt) + (D,N,Qp): D' = Dt D, N' = D Nt + Dt N, Q' = D^2 Qt + Dt^2 Qp
                Nn = [core.add(core.mul(D, Nt[i]), core.mul(Dt, N[i])) for i in range(6)]
                Qn = core.add(core.mul(core.mul(D, D), Qt), core.mul(core.mul(Dt, Dt), Qp))
                Dt, Nt, Qt = core.mul(Dt, D), Nn, Qn
            ok, info = core.verify_witness(F_, A_, 1, Dt, Nt, Qt)
            print(f'COMBINED: d={core.degree(Dt)} verify={ok} {info}', flush=True)
            with open(f'{WD}/bundle-pieces.json', 'w') as f:
                json.dump([{'D': core.to_str(D), 'A_S': core.to_str(AS), 'N': [core.to_str(n) for n in N],
                            'Q': core.to_str(Qp)} for (N, Qp, AS), D in zip(pieces, Ds)], f, indent=1)
            if ok:
                qcp.save_witness(Dt, Nt, Qt)
                print('WITNESS SAVED to witness-candidate.json', flush=True)
            sys.exit(0)
        else:
            print('   rational reconstruction verified false; continuing with more primes', flush=True)
print('lift failed after primes', flush=True)
