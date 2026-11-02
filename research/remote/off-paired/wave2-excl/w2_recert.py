"""Multi-prime kernel-lift certification for receipts marked kernel-lift-incomplete: compute the rref'd modular kernel basis at
several primes (same pivot structure), CRT the entries, rationally reconstruct, verify EXACTLY.
usage: python3 w2_recert.py --labels 'X1*M*H,X1*S*FY0' --job NAME [--primes p1,p2,p3,p4]"""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from w2lib import char_of
from hodge import core, symmetry, fixed_d
from hodge.core import NV, mul, sub, power
from fractions import Fraction

def kernel_basis(F, A, q, D, p, G):
    S, sub_ = symmetry.restricted_system(F, A, q, D, p, G)
    C, R = S.ncols, S.nrows
    compressed = R * (C + 1) > 40_000_000
    Maug, k = w2lib._compressed_matrix(S) if compressed else (w2lib._build_modular(S), R)
    Mr, rank_aug = Maug.rref()
    piv = []
    for i in range(rank_aug):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        piv.append(j)
    consistent = not (piv and piv[-1] == C)
    vecs = w2lib.kernel_vectors(S, Maug)
    return S, vecs, rank_aug - (0 if consistent else 1), rank_aug, consistent, compressed

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--labels', required=True); ap.add_argument('--job', required=True)
    ap.add_argument('--primes', default=','.join(str(p) for p in core.PRIMES[1:5]))
    a = ap.parse_args()
    primes = [int(x) for x in a.primes.split(',')]
    F, A, q = w2lib.load(); G = w2lib.group(F)
    for label in a.labels.split(','):
        t0 = time.time()
        D = w2lib.parse_expr(label)
        bases = []; info = None
        for p in primes:
            S, vecs, rM, rA, cons, comp = kernel_basis(F, A, q, D, p, G)
            if info is None: info = (S, rM, rA, cons, comp)
            bases.append((p, vecs))
            print(f'[{label}] p={p} rank_M={rM} rank_aug={rA} consistent={cons} kernel {len(vecs)}', flush=True)
        S, rM, rA, cons, comp = info
        C = S.ncols
        n = len(bases[0][1])
        if any(len(v) != n for _, v in bases):
            rec = {'exp': 'fixedD-recert', 'label': label, 'D': label, 'd': core.degree(D), 'char': char_of(D), 'job': a.job, 'primes': primes,
                   'result': 'kernel-dims-differ-between-primes', 'dims': [len(v) for _, v in bases], 'seconds': round(time.time() - t0, 1)}
            w2lib.log(rec); print(rec, flush=True); continue
        # CRT entrywise
        Fq1 = power(F, q + 1); AD2 = mul(A, mul(D, D))
        nver = 0; bad = 0; maxh = 0
        for i in range(n):
            rat = []; ok = True
            for j in range(C):
                acc, m = bases[0][1][i][j], bases[0][0]
                for p, vecs in bases[1:]:
                    acc, m = core.crt_pair(acc, m, vecs[i][j], p)[0], m * p
                r = core.ratrecon(acc, m) if acc else Fraction(0)
                if r is None: ok = False; break
                rat.append(r)
            if not ok: bad += 1; continue
            N = [{} for _ in range(NV)]; Q = {}
            for (kind, ii, mm), x in zip(S.cols, rat):
                if x == 0: continue
                if kind == 'N': N[ii][mm] = x
                else: Q[mm] = x
                maxh = max(maxh, abs(x.numerator), x.denominator)
            E = sub(sub(core.defect(F, A, D, q, N), AD2), mul(Fq1, Q))
            if not E: nver += 1
        proven = (not cons) and nver == n == C - rM
        rec = {'exp': 'fixedD-recert', 'label': label, 'D': label, 'd': core.degree(D), 'char': char_of(D), 'job': a.job, 'primes': primes,
               'cols': C, 'rows': S.nrows, 'rank_M': rM, 'rank_aug': rA, 'nullity': C - rM, 'kernel_verified': nver, 'kernel_total': n,
               'kernel_unrecon': bad, 'kernel_max_height': maxh, 'consistent': cons, 'compressed': comp,
               'result': 'PROVEN-inconsistent' if proven else ('solvable-modp' if cons else 'inconsistent-modp(kernel-lift-incomplete)'),
               'rigor': 'kernel-certified(multi-prime)' if proven else None, 'seconds': round(time.time() - t0, 1)}
        w2lib.log(rec); print(json.dumps(rec, default=str), flush=True)
        with open(f'{w2lib.PK}/receipts_{a.job}.jsonl', 'a') as f: f.write(json.dumps(rec, default=str) + '\n')

if __name__ == '__main__':
    main()
