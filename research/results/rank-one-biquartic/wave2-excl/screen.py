"""Normalized-E2 screening (mod p, full Koszul correction space + full C space) of E1-admissible denominators.
Usage:
  python3 screen.py gens d parity prime            -- every non-Jacobian kernel direction (basis of K_d / J_d, F-reduced family),
                                                      plus the Jacobian-part basis sample
  python3 screen.py random d parity prime n seed   -- n random members of the full F-reduced E1 family
  python3 screen.py randomNJ d parity prime n seed -- n random members of (non-Jacobian span) only
  python3 screen.py expr prime 'D expr' label      -- a specific rational D (E1 lift computed first)
Each experiment is logged to wave2-excl.log.jsonl."""
import sys, json, time, random
sys.path.insert(0, '.')
import w2
from flint import nmod_mat

mode = sys.argv[1]


def family_split(R, Fp, Ap, Fi, d, pD):
    """F-reduced E1 family; split basis into Jacobian part and a complement (non-Jacobian directions mod J)."""
    fam = w2.e1_family(R, Fp, Ap, Fi, d, pD, dpred=w2.not_div_X06)
    basis = fam['basis']
    if not basis:
        return fam, [], []
    p = R.p
    idx = w2.Index(d, pD, w2.not_div_X06)
    # Jacobian generators restricted to F-reduced monomials: J_d cap W_d = image of dF(G) projected? Not a subspace of W_d in general.
    # Instead: J_d cap W_d = {dF(G) : dF(G) has no X0^6-divisible monomial}. Compute via nullspace of the X0^6-divisible rows.
    gens = []
    for i in range(w2.NV):
        pi = None if pD is None else (pD - (0 if i < 3 else 1)) % 2
        for e in w2.monomials(d - 5, pi):
            gens.append(Fi[i] * R.mono(e))
    allidx = w2.Index(d, pD)
    bad = [j for j, e in enumerate(allidx.mons) if e[0] >= 6]
    n = len(gens)
    flat = [0] * (len(bad) * n)
    badpos = {j: r for r, j in enumerate(bad)}
    for c, P in enumerate(gens):
        for e, v in P.to_dict().items():
            j = allidx.idx[e]
            if j in badpos:
                flat[badpos[j] * n + c] = v
    Mb = nmod_mat(len(bad), n, flat, p)
    X, nul = Mb.nullspace()
    Jred = []
    for k in range(nul):
        P = R.zero()
        for c in range(n):
            v = int(X[c, k])
            if v:
                P += gens[c] * v
        if not P.is_zero():
            Jred.append(P)
    # rank of Jred and complement of span(Jred) inside span(basis)
    k = len(basis)
    flatJ = [0] * (len(Jred) * idx.n)
    for r_, P in enumerate(Jred):
        for e, v in P.to_dict().items():
            flatJ[r_ * idx.n + idx.idx[e]] = v
    MJ = nmod_mat(len(Jred), idx.n, flatJ, p) if Jred else None
    rJ = MJ.rank() if MJ is not None else 0
    # greedy complement: add basis vectors that increase the rank
    cur = [list(row) for row in (MJ.tolist() if MJ is not None else [])]
    nonJ = []
    currank = rJ
    for (D, nn) in basis:
        row = idx.vec(D)
        test = nmod_mat(len(cur) + 1, idx.n, [x for rr in cur for x in rr] + row, p)
        rk = test.rank()
        if rk > currank:
            nonJ.append((D, nn)); cur.append(row); currank = rk
        if currank == k:
            break
    return fam, Jred, nonJ


def run(R, Fp, Fi, B, D, n, label, d, pD, p, note=''):
    t0 = time.time()
    r = B.solve(D, n)
    rec = dict(exp='e2-screen', p=p, d=d, parity=pD, D=label, consistent=r['consistent'], rank_M=r['rank_M'],
               rank_aug=r['rank_aug'], rows=r['rows'], cols=r['cols'], n_koszul=r['n_koszul'], n_C=r['n_C'],
               seconds=round(time.time() - t0, 1), note=note)
    if r['consistent']:
        K, Cp = B.unpack(r['x'])
        N = [n[i] + K[i] for i in range(w2.NV)]
        rec['WITNESS_CANDIDATE_modp'] = dict(D=str(D), N=[str(x) for x in N], Q=str(Cp))
        print('!!!!!!!! CONSISTENT (mod p) !!!!!!!!', label, file=sys.stderr, flush=True)
        json.dump(rec, open(f'WITNESS_CANDIDATE_{int(time.time())}.json', 'w'), default=str)
    w2.log(rec)
    print(json.dumps({k: v for k, v in rec.items() if k != 'WITNESS_CANDIDATE_modp'}), flush=True)
    return r


if mode in ('gens', 'random', 'randomNJ'):
    d, par, p = int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    pD = None if par == 'm' else int(par)
    R, Fp, Ap, Fi, F, A = w2.load(p)
    fam, Jred, nonJ = family_split(R, Fp, Ap, Fi, d, pD)
    print(f'd={d} par={par} p={p}: F-reduced E1 family dim {fam["dimK"]}, Jacobian part dim {len(Jred)}, non-Jacobian directions {len(nonJ)}', flush=True)
    w2.log(dict(exp='e1-family-split', p=p, d=d, parity=par, dimK=fam['dimK'], dimJ_reduced=len(Jred), n_nonJ=len(nonJ)))
    B = w2.E2Builder(R, Fp, Fi, d, pD)
    basis = fam['basis']
    if mode == 'gens':
        for k, (D, n) in enumerate(nonJ):
            run(R, Fp, Fi, B, D, n, f'nonJ-direction {k} of K_{d}^{par} (mod p rep, F-reduced)', d, pD, p,
                note='representative of a non-Jacobian annihilator direction (basis element of K cap W mod J-part)')
        # also test a few random members of the Jacobian part alone
        random.seed(1)
        idx = w2.Index(d, pD, w2.not_div_X06)
        for k in range(3):
            D = R.zero()
            for P in Jred:
                D += P * random.randint(1, p - 1)
            n = w2.e1_lift(R, Fp, Ap, Fi, D)
            run(R, Fp, Fi, B, D, n, f'random member #{k} of the Jacobian part J_{d}^{par} cap W_{d} (mod p)', d, pD, p,
                note='pure polar denominator dF(G), F-reduced')
    else:
        nn, seed = int(sys.argv[5]), int(sys.argv[6])
        random.seed(seed)
        src = basis if mode == 'random' else nonJ
        for k in range(nn):
            D = R.zero(); n = [R.zero() for _ in range(w2.NV)]
            for (Dk, nk) in src:
                c = random.randint(1, p - 1)
                D += Dk * c
                for i in range(w2.NV):
                    n[i] += nk[i] * c
            run(R, Fp, Fi, B, D, n, f'random member #{k} (seed {seed}) of {"full F-reduced E1 family" if mode == "random" else "non-Jacobian span"} d={d} par={par}', d, pD, p)
elif mode == 'expr':
    p = int(sys.argv[2]); expr = sys.argv[3]; label = sys.argv[4]
    R, Fp, Ap, Fi, F, A = w2.load(p)
    from hodge import core
    D = R.from_dict(core.parse(expr))
    d = D.total_degree(); pD = R.parity(D)
    n = w2.e1_lift(R, Fp, Ap, Fi, D)
    if n is None:
        rec = dict(exp='e2-screen', p=p, d=d, parity=pD, D=label, e1=False, consistent=False, note='E1 fails: A D not in J(F) mod p')
        w2.log(rec); print(json.dumps(rec)); sys.exit(0)
    B = w2.E2Builder(R, Fp, Fi, d, pD)
    run(R, Fp, Fi, B, D, n, label, d, pD, p, note='explicit denominator; E1 lift found')
