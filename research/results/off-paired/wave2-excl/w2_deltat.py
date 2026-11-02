"""Full fixed-D solve of the cubic-cover discriminant Delta_t = 4 Qx^3 + 27 f6^2 g6 (degree 18, G-invariant, irreducible),
chunked row compression (memory-light), rank sandwich + kernel-lift certification.  usage: python3 w2_deltat.py JOB [label expr]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from hodge import core, fixed_d, symmetry
from hodge.core import mul, add, sub, scale, power
import numpy as np
from flint import nmod_mat

job = sys.argv[1]
label = sys.argv[2] if len(sys.argv) > 2 else 'Delta_t'
expr = sys.argv[3] if len(sys.argv) > 3 else '4*Qx**3+27*f6**2*g6'
F, A, q = w2lib.load(); G = w2lib.group(F)
D = w2lib.parse_expr(expr)
print(f'[{job}] D={label} degree {core.degree(D)} terms {len(D)} char {w2lib.char_of(D)}', flush=True)
e1 = w2lib.e1_receipt(F, A, q, D, w2lib.P1, label, 'discriminant'); e1['job'] = job
w2lib.log(e1); print(e1, flush=True)
if not e1['e1']:
    sys.exit(0)
for p in (w2lib.P1, w2lib.P2):
    t0 = time.time()
    S, sub = symmetry.restricted_system(F, A, q, D, p, G, verbose=True)
    M = fixed_d._sparse_matrix(S)
    C = S.ncols; k = min(S.nrows, C + 1 + 40)
    MT = M.T.tocsr(); del M
    rng = np.random.default_rng(1)
    out = np.empty((k, C + 1), dtype=np.int64)
    blk = 512
    for s in range(0, k, blk):
        e = min(k, s + blk)
        R = rng.integers(0, p, size=(e - s, S.nrows), dtype=np.int64)
        lo = (MT @ (R & 0xFFFF).T) % p
        hi = (MT @ (R >> 16).T) % p
        out[s:e] = ((lo + (hi * (1 << 16)) % p) % p).T
        del R, lo, hi
    t1 = time.time(); print(f'[{job}] compressed k={k} C={C} in {t1-t0:.0f}s', flush=True)
    del MT
    Maug = nmod_mat(k, C + 1, out.ravel().tolist(), p); del out
    Mr, rank_aug = Maug.rref()
    t2 = time.time()
    pivots = []
    for i in range(rank_aug):
        j = 0
        while j <= C and int(Mr[i, j]) == 0: j += 1
        pivots.append(j)
    consistent = not (pivots and pivots[-1] == C)
    rank_M = rank_aug - (0 if consistent else 1)
    e = w2lib.explicit_kernel_dim(F, A, q, D, G)
    rec = {'exp': 'fixedD', 'D': label, 'family': 'discriminant', 'd': core.degree(D), 'char': w2lib.char_of(D), 'prime': p,
           'rows': S.nrows, 'cols': C, 'nN': S.nN, 'nQ': S.nQ, 'k_rows': k, 'rank_M': rank_M, 'rank_aug': rank_aug,
           'nullity': C - rank_M, 'explicit_kernel': e, 'consistent': consistent, 'compressed': True, 'job': job,
           't_rref': round(t2 - t1, 1)}
    print(f'[{job}] p={p} consistent={consistent} rank_M={rank_M} rank_aug={rank_aug} nullity={C-rank_M} explicit={e} rref {t2-t1:.0f}s', flush=True)
    if consistent:
        rec['result'] = 'solvable-modp'
        x = [0] * C
        for i, j in enumerate(pivots): x[j] = int(Mr[i, C])
        N, Q = S.unpack(x)
        json.dump({'D': core.to_str(D), 'p': p, 'N': [{str(k_): v for k_, v in n.items()} for n in N], 'Q': {str(k_): v for k_, v in Q.items()}},
                  open(f'/tmp/hodge/results/off-paired/wave2-excl/{job}_solution_modp.json', 'w'))
        print('!!!!! CONSISTENT mod p -- attempting exact lift', flush=True)
        st, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, primes=[p] + [pp for pp in core.PRIMES if pp != p], verbose=True, compressed=True)
        print('lift status', st, flush=True)
        if st == 'solved':
            json.dump({'instance': 'off-paired', 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])},
                      open('/tmp/hodge/results/off-paired/witness-candidate.json', 'w'))
            rec['result'] = 'WITNESS'
    elif rank_M == C - e:
        rec['result'] = 'PROVEN-inconsistent'; rec['rigor'] = 'explicit-kernel'
    else:
        vecs = w2lib.kernel_vectors(S, Maug)
        del Maug, Mr
        nver, ntot, maxh, bad = w2lib.certify_kernel(F, A, q, D, S, vecs, p, verbose=True)
        rec.update({'kernel_verified': nver, 'kernel_total': ntot, 'kernel_max_height': maxh, 'kernel_unrecon': bad})
        rec['result'] = 'PROVEN-inconsistent' if nver == ntot == C - rank_M else 'inconsistent-modp(kernel-lift-incomplete)'
        if nver == ntot == C - rank_M: rec['rigor'] = 'kernel-certified'
    rec['seconds'] = round(time.time() - t0, 1)
    w2lib.log(rec); print(rec, flush=True)
    with open(f'/tmp/hodge/results/off-paired/wave2-excl/receipts_{job}.jsonl', 'a') as f:
        f.write(json.dumps(rec, default=str) + '\n')
    if consistent:
        break
