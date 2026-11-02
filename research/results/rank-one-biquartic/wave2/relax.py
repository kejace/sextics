"""Jacobian-quotient relaxation of E2 on an E1-admissible family (exact mod p).

For D(a) = sum a_m D_m with lifts N0(a) = sum a_m N_m, E2 requires B_{D(a)}(N0(a)) in J(F)_{2d} (necessary;
"Hessian relaxation": for Jacobian-type D = dF(w) it reads A w^T Hess(F) w in J).  The map
  Phi: Sym^2(family) -> (S/J)_{2d}^{pA},  a (x) a -> [B_D(N0)]
is linear in the pair coordinates; if rank Phi = #pairs there is no nonzero solution over F_p-bar, and by rank
semicontinuity (J-ranks equal the regular-sequence values) none over C.

  python3 wave2/relax.py --d 6 --pD 0 [--msolve-out file.ms] [--tag relax]
"""
import sys, os, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
from flint import nmod_mat

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--pD', type=int, default=0)
ap.add_argument('--prime-index', type=int, default=0)
ap.add_argument('--tag', default='relax')
ap.add_argument('--msolve-out', default=None)
ap.add_argument('--kernel-json', default=None, help='out/e1kernel_d{d}_p{pD}.json (kernel reps + lifts mod p)')
args = ap.parse_args()

F, A, q = load()
d, pD = args.d, args.pD
p = core.PRIMES[args.prime_index]
Ap = core.reduce_mod(A, p); Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
pA = core.parity_of(A)
t0 = time.time()
# family: J-type generators (D = m F_i, N = A m e_i) + kernel reps with lifts
fam = []
for i in range(NV):
    par = (pD - par_idx(i)) % 2
    for m in monomials(d - 5, par):
        Dm = shift(Fi[i], m, 1, p)
        Nm = [{} for _ in range(NV)]; Nm[i] = shift(Ap, m, 1, p)
        fam.append(('J', i, m, Dm, Nm))
nJ = len(fam)
kpath = args.kernel_json or f'{OUT}/e1kernel_d{d}_p{pD}.json'
nker = 0
if os.path.exists(kpath):
    kj = json.load(open(kpath))
    assert kj['p'] == p, 'kernel file prime mismatch'
    for b, lift in zip(kj['basis_modp'], kj['lifts_modp']):
        Dm = core.parse(b, p)
        Nm = [core.parse(s, p) if s != '0' else {} for s in lift]
        fam.append(('K', None, None, Dm, Nm)); nker += 1
else:
    print(f'[relax] no kernel file {kpath}; using J-type family only', file=sys.stderr)
Kf = len(fam)
# independence of the D's (mod p) -- for d <= 9 the J-type generators are independent
SD = Space(d, pD)
flat = [0] * (Kf * SD.n)
for r_, (_, _, _, Dm, _) in enumerate(fam):
    for e, v in Dm.items():
        flat[r_ * SD.n + SD.idx[e]] = v
rankD = nmod_mat(Kf, SD.n, flat, p).rank()
print(f'[relax] d={d} pD={pD} family: {nJ} J-type + {nker} kernel = {Kf}; rank of D-span = {rankD}', file=sys.stderr, flush=True)
# verify lifts
for (_, _, _, Dm, Nm) in fam:
    NF = {}
    for i in range(NV):
        NF = add(NF, mul(Nm[i], Fi[i], p), p)
    assert sub(NF, mul(Ap, Dm, p), p) == {}, 'lift check failed'
# reducer in degree 2d, parity pA
t1 = time.time()
Srow = Space(2 * d, pA)
Jm, mons = jacobian.jacobian_matrix(F, 2 * d, p, pA, Srow.mons)
Jr, rankJ = Jm.rref()
ncol = Srow.n
pivots = []
for i in range(rankJ):
    j = 0
    while j < ncol and int(Jr[i, j]) == 0:
        j += 1
    pivots.append(j)
pivset = set(pivots)
std = [j for j in range(ncol) if j not in pivset]
std_pos = {j: t for t, j in enumerate(std)}
ring_dim = len(std)
# R_std: rank x n_std block of rref rows restricted to standard columns
Rstd = nmod_mat(rankJ, ring_dim, [int(Jr[i, j]) for i in range(rankJ) for j in std], p)
print(f'[relax] J_{2*d}^{pA}: gens {Jm.nrows()} rank {rankJ} ring_dim {ring_dim} ({time.time()-t1:.0f}s)', file=sys.stderr, flush=True)
# pair products
t2 = time.time()
pairs = [(m, n) for m in range(Kf) for n in range(m, Kf)]
npairs = len(pairs)
# rows = pairs, columns = Srow monomials; then NF = V_std - V_piv Rstd
Vflat_piv = [0] * (npairs * rankJ)
Vflat_std = [0] * (npairs * ring_dim)
piv_pos = {j: t for t, j in enumerate(pivots)}
for t_, (m, n) in enumerate(pairs):
    Bmn = B_of(fam[m][3], fam[n][4], p)
    if m != n:
        Bmn = add(Bmn, B_of(fam[n][3], fam[m][4], p), p)
    for e, v in Bmn.items():
        j = Srow.idx[e]
        if j in pivset:
            Vflat_piv[t_ * rankJ + piv_pos[j]] = v
        else:
            Vflat_std[t_ * ring_dim + std_pos[j]] = v
    if t_ % 500 == 0:
        print(f'  pair {t_}/{npairs} ({time.time()-t2:.0f}s)', file=sys.stderr, flush=True)
Vp = nmod_mat(npairs, rankJ, Vflat_piv, p)
Vs = nmod_mat(npairs, ring_dim, Vflat_std, p)
NF = Vs - Vp * Rstd
rankPhi = NF.rank()
info = {'batch': 'relax', 'd': d, 'pD': pD, 'p': p, 'family_J': nJ, 'family_kernel': nker, 'K': Kf, 'rankD': int(rankD),
        'npairs': npairs, 'ring_dim_2d': ring_dim, 'rank_Phi': int(rankPhi),
        'excluded_over_Fp': bool(rankPhi == npairs), 'time_s': round(time.time() - t0, 1)}
print(json.dumps(info), flush=True)
log_line(f'{OUT}/{args.tag}.jsonl', info)
if args.msolve_out:
    # quadrics: for each standard monomial row (after transposing), the quadric sum_{pairs} NF[pair, row] a_m a_n
    # write only a basis of the row space (rank many) to keep the file small
    NFt = NF.transpose()               # ring_dim x npairs
    Rr, rk = NFt.rref()
    names = [f'a{m}' for m in range(Kf)]
    with open(args.msolve_out, 'w') as fh:
        fh.write(','.join(names) + '\n' + f'{p}\n')
        polys = []
        for i in range(rk):
            terms = []
            for t_, (m, n) in enumerate(pairs):
                v = int(Rr[i, t_])
                if v:
                    terms.append(f'{v}*{names[m]}*{names[n]}' if m != n else f'{v}*{names[m]}^2')
            polys.append('+'.join(terms))
        fh.write(',\n'.join(polys) + '\n')
    print(f'[relax] wrote {rk} quadrics in {Kf} variables to {args.msolve_out}', file=sys.stderr)
