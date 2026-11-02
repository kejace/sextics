"""Step 4a: E1-admissible 'cone' subfamilies.  For each degree d, variable subset V (|V| = 2,3,4: cones with
vertex a plane/line over curves/surfaces), and character chi, compute dim( R_d^chi[V] ∩ (J:A)_d ) exactly
(mod p rank test on the rational E1 basis), excluding the Euler gauge.  Writes wave2-cones.jsonl; nonzero
subfamilies are the inputs of w2_cones_search.py."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/off-paired')
import w2lib as W
from fractions import Fraction
import flint
F, A, q = W.load()
d0, d1 = int(sys.argv[1]), int(sys.argv[2])
out = open('/tmp/hodge/results/off-paired/wave2-cones.jsonl', 'a')
names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
t00 = time.time(); nfound = 0; ntest = 0
for d in range(d0, d1 + 1):
    for par in (0, 1):
        for c in W.all_chars_parity(par):
            S = W.SectorSpaces(F, A, q, d, c)
            if S.dimE1 - S.dim_gauge <= 0:
                continue
            Dm = S.Dmons; K = len(Dm); idx = {m: j for j, m in enumerate(Dm)}
            # E1 basis rows mod p, gauge rows mod p
            e1rows = [S._modrow([Fraction(v) for v in row]) for row in S.E1]
            grows = [S._modrow(S._vec(P, Dm)) for P in S.gauge]
            rE1 = S._rank_rows(e1rows, K)
            for size in (2, 3, 4):
                for V in itertools.combinations(range(6), size):
                    cols = [j for j, m in enumerate(Dm) if all(m[i] == 0 for i in range(6) if i not in V)]
                    if not cols:
                        continue
                    ntest += 1
                    vrows = []
                    for j in cols:
                        r = [0] * K; r[j] = 1; vrows.append(r)
                    rV = len(cols)
                    rsum = S._rank_rows(e1rows + vrows, K)
                    dim_int = rE1 + rV - rsum
                    if dim_int <= 0:
                        continue
                    # remove the part inside the gauge: dim(V ∩ E1 ∩ gauge) = dim(V∩gauge) since gauge ⊂ E1
                    rG = S._rank_rows(grows, K)
                    rsumG = S._rank_rows(grows + vrows, K) if grows else rV
                    dim_gauge_int = rG + rV - rsumG
                    net = dim_int - dim_gauge_int
                    if net <= 0:
                        continue
                    # is the intersection inside J (Jacobian type)?  dim(V ∩ (gauge+J))
                    jrows = grows + [S._modrow(S._vec(P, Dm)) for P in S.jac_basis]
                    rJ = S._rank_rows(jrows, K)
                    rsumJ = S._rank_rows(jrows + vrows, K)
                    dim_J_int = rJ + rV - rsumJ
                    rec = {'d': d, 'chi': list(c), 'vars': [names[i] for i in V], 'dim_V': rV, 'dim_cap_E1': dim_int,
                           'dim_cap_gauge': dim_gauge_int, 'net': net, 'dim_cap_jac_plus_gauge': dim_J_int,
                           'non_jacobian': dim_int - dim_J_int}
                    out.write(json.dumps(rec) + '\n'); out.flush(); nfound += 1
                    print(json.dumps(rec), flush=True)
    W.log_batch({'step': 'cones-dims', 'd': d, 'subsets_tested_cum': ntest, 'nonzero_cum': nfound,
                 'seconds_cum': round(time.time() - t00, 1), 'detail_file': 'wave2-cones.jsonl'})
print('done', nfound, 'nonzero subfamilies of', ntest, 'tests in', round(time.time() - t00), 's')
