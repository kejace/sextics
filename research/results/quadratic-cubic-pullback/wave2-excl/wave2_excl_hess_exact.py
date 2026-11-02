"""Exact completion of the Hessian relaxation: lift the (small-coefficient) mod-p kernel basis of
      W -> [A sum_{a<=b} W_ab m_a m_b F_{i_a i_b}]  in  S_{2d}/J(F)_{2d}
to Q, certify each lifted vector exactly (membership in J(F)_{2d} by exact linear algebra: Dixon solve on a
maximal-rank square subsystem + full exact verification), then compute the exact annihilator quadrics and the
candidate locus {w : w w^T in K_H} with Singular over Q (dimension + minimal associated primes).

usage: python3 wave2_excl_hess_exact.py <jobname> <d> <parity> <prime of the KH file>
"""
import sys, os, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/quadratic-cubic-pullback')
import w2x
from w2x import core, F, A, monomials, mul, add, sub, scale, shift, pderiv, logline, dim_J
from fractions import Fraction
from math import gcd
from flint import fmpz_mat, fmpq_mat, fmpq, nmod_mat
JOB, D_DEG, PAR, P = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
LOG = os.path.join(w2x.LANE, f'wave2-excl.{JOB}.jsonl')
parz = [0, 0, 0, 1, 1, 1]
pA = core.parity_of(A)


def main():
    t0 = time.time()
    d = D_DEG
    K = json.load(open(os.path.join(w2x.LANE, f'wave2-excl-hess-hessL{d}-d{d}-par{PAR}-p{P}-KH.json')))
    wb = [(i, tuple(m)) for i, m in K['wb']]
    Wkeys = [tuple(x) for x in K['Wkeys']]
    basis = [{int(a): v for a, v in row.items()} for row in K['basis']]
    nw, nW, rk = len(wb), len(Wkeys), len(basis)
    Fi = [pderiv(F, i) for i in range(6)]
    Fij = [[pderiv(Fi[i], j) for j in range(6)] for i in range(6)]
    # exact target polynomial for each lifted basis vector
    def target(row):
        T = {}
        for a, v in row.items():
            (ia, ma), (ib, mb) = wb[Wkeys[a][0]], wb[Wkeys[a][1]]
            t = mul(A, shift(Fij[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 1))
            if Wkeys[a][0] != Wkeys[a][1]:
                t = scale(t, 2)
            T = add(T, scale(t, v))
        return T
    targets = [target(row) for row in basis]
    # split by Y-parity; verify each parity part in its own Jacobian block (J is tau-stable)
    def part(T, par):
        return {e: v for e, v in T.items() if (e[3] + e[4] + e[5]) % 2 == par}
    ok_all = True
    rank_info = {}
    p = core.PRIMES[1]
    t1 = time.time()
    for tpar in (0, 1):
        tg = [part(T, tpar) for T in targets]
        if not any(tg):
            continue
        gens = []
        for i in range(6):
            for m in monomials(2 * d - 5, (tpar - parz[i]) % 2):
                gens.append(shift(Fi[i], m, 1))
        mons = monomials(2 * d, tpar)
        col = {e: j for j, e in enumerate(mons)}
        ng, nm = len(gens), len(mons)
        flat = [0] * (ng * nm)
        for r, g in enumerate(gens):
            for e, v in g.items():
                flat[r * nm + col[e]] = (v.numerator * pow(v.denominator, -1, p)) % p
        M = nmod_mat(ng, nm, flat, p)
        Mr, rank = M.rref()
        assert rank == dim_J(2 * d, tpar), ('J rank not Hilbert', rank, dim_J(2 * d, tpar))
        pivcols = []
        for i in range(rank):
            j = 0
            while int(Mr[i, j]) == 0:
                j += 1
            pivcols.append(j)
        Mt = M.transpose(); Mtr, rank2 = Mt.rref()
        pivrows = []
        for i in range(rank2):
            j = 0
            while int(Mtr[i, j]) == 0:
                j += 1
            pivrows.append(j)
        del M, Mr, Mt, Mtr
        pc_index = {c: ii for ii, c in enumerate(pivcols)}
        S = fmpz_mat(rank, rank); scal = []
        for jj, r in enumerate(pivrows):
            g = gens[r]; den = 1
            for v in g.values():
                den = den * v.denominator // gcd(den, v.denominator)
            scal.append(den)
            for e, v in g.items():
                c = col[e]
                if c in pc_index:
                    S[pc_index[c], jj] = int(v * den)
        # check: square subsystem must be nonsingular mod p (otherwise choose differently)
        B = fmpz_mat(rank, rk); tden = []
        for k, T in enumerate(tg):
            den = 1
            for v in T.values():
                den = den * v.denominator // gcd(den, v.denominator)
            tden.append(den)
            for e, v in T.items():
                c = col[e]
                if c in pc_index:
                    B[pc_index[c], k] = int(v * den)
        Xq = S.solve(B)
        nbad = 0
        for k in range(rk):
            if not tg[k]:
                continue
            acc = {}
            for jj, r in enumerate(pivrows):
                xv = Fraction(int(Xq[jj, k].p), int(Xq[jj, k].q))
                if xv == 0:
                    continue
                acc = add(acc, scale(gens[r], xv * scal[jj]))
            resid = {e: v for e, v in sub(acc, scale(tg[k], tden[k])).items() if v != 0}
            if resid:
                nbad += 1
        rank_info[tpar] = dict(J_rank=rank, n_targets=sum(1 for T in tg if T), failures=nbad)
        if nbad:
            ok_all = False
    t2 = time.time()
    rank = rank_info
    rec = dict(exp='hess-exact', d=d, parity=PAR, prime_of_basis=P, dim_KH=rk, nw=nw, nW=nW, J_blocks=rank_info,
               exact_membership_all=ok_all, seconds_solve=round(t2 - t1, 1), job=JOB)
    print(json.dumps(rec), flush=True)
    if not ok_all:
        logline(LOG, rec); return
    # exact annihilator quadrics
    Bz = fmpz_mat(rk, nW, [basis[i].get(a, 0) for i in range(rk) for a in range(nW)])
    Lam, nlam = Bz.nullspace()
    quadrics = []
    for k in range(nlam):
        colv = [int(Lam[idx, k]) for idx in range(nW)]
        g = 0
        for c in colv:
            g = gcd(g, abs(c))
        g = g or 1
        terms = []
        for idx, (a, b) in enumerate(Wkeys):
            c = colv[idx] // g
            if c:
                terms.append(f'{c}*w{a}*w{b}')
        quadrics.append('+'.join(terms).replace('+-', '-'))
    rec['n_quadrics_exact'] = nlam
    # Singular over Q: dimension of the cone and minimal associated primes
    sing = os.path.join(w2x.LANE, f'wave2-excl-hess-exact-{JOB}-d{d}-par{PAR}.sing')
    with open(sing, 'w') as f:
        f.write('ring r = 0, (' + ','.join(f'w{a}' for a in range(nw)) + '), dp;\n')
        f.write('ideal I = ' + ',\n'.join(quadrics) + ';\n')
        f.write('ideal G = std(I);\nprint("DIM " + string(dim(G)));\nLIB "primdec.lib";\nlist L = minAssGTZ(I);\nprint("NCOMP " + string(size(L)));\n')
        f.write('int k; for (k = 1; k <= size(L); k++) { print("COMPONENT " + string(k) + " dim " + string(dim(std(L[k])))); print(L[k]); }\nquit;\n')
    try:
        r = subprocess.run(['Singular', '-q', sing], capture_output=True, text=True, timeout=600)
        out = r.stdout
    except subprocess.TimeoutExpired:
        out = 'TIMEOUT'
    rec['singular_output'] = out[-3000:]
    rec['seconds'] = round(time.time() - t0, 1)
    logline(LOG, rec)
    print(out[-2500:])
    print(json.dumps({k: v for k, v in rec.items() if k != 'singular_output'}), flush=True)


if __name__ == '__main__':
    main()
