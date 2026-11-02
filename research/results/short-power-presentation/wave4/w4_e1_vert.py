"""E1 kernels of structured denominator families (vertical ring, cones) by the T-reduction, two primes.
Usage: python3 w4_e1_vert.py <out.jsonl> <family: vert|coneXM|coneLY> d1 d2 ..."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/wave4')
from w4_tlib import *
from w2c_lib import fam_vert, fam_coneXM, fam_coneLY
FAMS = {'vert': fam_vert, 'coneXM': fam_coneXM, 'coneLY': fam_coneLY}
out_path, fname = sys.argv[1], sys.argv[2]
degs = [int(x) for x in sys.argv[3:]]
F, A, q, rec = core.load_instance(IID)
P1, P2 = core.PRIMES[0], core.PRIMES[1]

def complement(Jvecs, Evecs, names, p):
    """E1 vectors reduced modulo span(Jvecs); returns list of nonzero remainders (independent, echelon form)."""
    K = len(names); nidx = {n: i for i, n in enumerate(names)}
    rows = []
    def tovec(v): 
        r = [0] * K
        for n, c in v.items(): r[nidx[n]] = c % p
        return r
    base = [tovec(v) for v in Jvecs]
    out = []
    for v in Evecs:
        cand = base + out + [tovec(v)]
        M = nmod_mat(len(cand), K, [x for r in cand for x in r], p)
        if M.rank() > len(base) + len(out):
            # reduce v modulo base+out via rref of the stacked matrix (last row carries v)
            Mr, rk = M.rref()
            # the reduced row for v is the one in the rref that is not in the span of base+out; recompute simply:
            Mb = nmod_mat(len(base) + len(out) + 1, K, [x for r in (base + out) for x in r] + tovec(v), p)
            # take v's reduction: solve via rref of [base+out] then eliminate pivots from v
            B, rb = nmod_mat(len(base) + len(out), K, [x for r in (base + out) for x in r], p).rref() if base + out else (None, 0)
            w = tovec(v)
            for i in range(rb):
                j = 0
                while int(B[i, j]) == 0: j += 1
                if w[j]:
                    c = w[j]
                    for k in range(j, K): w[k] = (w[k] - c * int(B[i, k])) % p
            # normalise leading entry to 1
            j = 0
            while j < K and w[j] == 0: j += 1
            inv = pow(w[j], -1, p)
            w = [(x * inv) % p for x in w]
            out.append(w)
    return [{names[i]: x for i, x in enumerate(w) if x} for w in out]

for d in degs:
    for par in (0, 1):
        fam = FAMS[fname](d, par)
        names = [n for n, _ in fam]
        if not fam:
            continue
        t0 = time.time()
        res = {}
        for p in (P1, P2):
            vecs, info = e1_kernel(F, A, fam, d, par, p)
            jv, infoj = in_J_family(F, fam, d, par, p)
            comp = complement(jv, vecs, names, p)
            Fp = reduce_mod(F, p)
            fdiv = [bool(jacobian.divisible_by(Fp, vec_to_poly(v, fam, p), p)[0]) if d >= 6 else False for v in comp]
            res[p] = {'kernel_dim': info['kernel_dim'], 'J_part_dim': infoj['kernel_dim'], 'interesting_dim': len(comp),
                      'hilbert_ok': info['hilbert_ok'], 'hilbert_ok_degd': infoj['hilbert_ok'],
                      'rank_J': info['rank_J'], 'ring_dim': info['ring_dim'], 'ring_dim_pred': info['ring_dim_pred'],
                      'rows': info['rows'], 'cols': info['cols'], 'fam_dim': info['K'], 'fam_rank': info['K'] - (info['cols'] - info['rank_total'] - (info['CJ'] - info['rank_J'])) if False else None,
                      'interesting': [sym(v, p) for v in comp], 'interesting_F_divides': fdiv,
                      't_rref': info['t_rref'] + infoj['t_rref']}
        # cross-prime rational reconstruction of the interesting vectors (canonical echelon form => comparable)
        recon = []
        c1, c2 = res[P1]['interesting'], res[P2]['interesting']
        if len(c1) == len(c2):
            for v1, v2 in zip(c1, c2):
                keys = set(v1) | set(v2); r = {}; ok = True
                for k in keys:
                    x, mod = core.crt_pair(v1.get(k, 0) % P1, P1, v2.get(k, 0) % P2, P2)
                    rr = core.ratrecon(x, mod)
                    if rr is None: ok = False; break
                    if rr: r[k] = str(rr)
                recon.append(r if ok else 'ratrecon-failed')
        rec_out = {'what': 'E1-kernel', 'family': fname, 'd': d, 'parD': par, 'fam_dim': len(fam),
                   'p1': P1, 'p2': P2, 'res_p1': res[P1], 'res_p2': res[P2], 'recon': recon,
                   'seconds': round(time.time() - t0, 1), 'file': out_path}
        with open(out_path, 'a') as f: f.write(json.dumps(rec_out) + '\n')
        log({'what': 'E1-kernel-' + fname, 'd': d, 'parD': par, 'prime': [P1, P2],
             'result': 'fam %d | kernel %d/%d | J-part %d/%d | interesting %d/%d | hilbert %s/%s' % (
                 len(fam), res[P1]['kernel_dim'], res[P2]['kernel_dim'], res[P1]['J_part_dim'], res[P2]['J_part_dim'],
                 res[P1]['interesting_dim'], res[P2]['interesting_dim'], res[P1]['hilbert_ok'], res[P2]['hilbert_ok']),
             'file': out_path})
        print('DONE', fname, d, par, 'fam', len(fam), 'kernel', res[P1]['kernel_dim'], res[P2]['kernel_dim'], 'J-part', res[P1]['J_part_dim'],
              'interesting', res[P1]['interesting_dim'], res[P2]['interesting_dim'], 'hilbert', res[P1]['hilbert_ok'], res[P2]['hilbert_ok'],
              'Fdiv', res[P1]['interesting_F_divides'], '%.0fs' % (time.time() - t0), flush=True)
print('ALL DONE', flush=True)
