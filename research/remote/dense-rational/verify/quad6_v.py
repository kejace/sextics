"""Independent degree-6 (and degree-5) Hessian-relaxation emptiness check.
d=6: D = gradF . w, w linear.  Sectors: even (w type +), odd (type -), mixed (both).  Condition w^T Hess(F) w in J_6
(valid since (J:A)_6 = J_6 in both parities).  Gauge w -> w + c z (type + only) fixed by dropping the X2-coefficient of w_{X2}.
Writes msolve input (prime p) with homogeneous quadrics; runs msolve -g 1; checks pure powers of every variable in the leading ideal.
d=5: exact rank over Q of Sym^2(C^6) -> S_4, S -> sum S_ij F_ij."""
import sys, time, json, subprocess
sys.path.insert(0, '/tmp/hodge/results/dense-rational/verify')
from vlib import *
p = int(sys.argv[1])
F, A, q = parse_rec()
Fp = reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(6)]
H = [[pderiv(Fi[i], j, p) for j in range(6)] for i in range(6)]
pz = [0, 0, 0, 1, 1, 1]
S, R = hilbert()
# --- degree 5 exact rank ---
FiQ = [pderiv(F, i) for i in range(6)]
HQ = [[pderiv(FiQ[i], j) for j in range(6)] for i in range(6)]
mons4 = monomials(4); idx4 = {m: k for k, m in enumerate(mons4)}
M5 = fmpq_mat(len(mons4), 21)
c = 0
for i in range(6):
    for j in range(i, 6):
        for e, v in HQ[i][j].items(): M5[idx4[e], c] = fmpq(v.numerator, v.denominator)
        c += 1
rk5 = M5.rank()
print('degree-5: rank of Sym^2 -> S_4 over Q:', rk5)
log('excl R-c: degree 5 excluded for every direction v (Sym^2(C^6) -> S_4, S -> sum S_ij F_ij injective, rank 21 exact)', 'CONFIRMED' if rk5 == 21 else 'REFUTED',
    f'exact fmpq rank {rk5} of 126 x 21', 'logic checked: E1 => N = A v + K gradF; E2 mod J => A v^T Hess v in J_10; (J:A)_4 = 0 both parities (sweep) => v^T Hess v = 0 => v v^T = 0 => v = 0; mixed v included since the quadratic form identity is checked componentwise in parity')
# --- degree 6 three sectors ---
results = {}
for sector in ('even', 'odd', 'mixed'):
    t0 = time.time()
    basis = []
    for i in range(6):
        for lm in monomials(1):
            typ = (pz[i] + yparity(lm)) % 2   # 0: type +, 1: type -
            if sector == 'even' and typ != 0: continue
            if sector == 'odd' and typ != 1: continue
            basis.append((i, lm))
    # gauge: w = c z  (type +): drop coordinate (2, X2)
    gauge_key = (2, (0, 0, 1, 0, 0, 0))
    keep = [b for b in basis if b != gauge_key]
    nvar = len(keep)
    # normal form mod J_6 (both parities): rref of J_6 generators
    rows6 = monomials(6); idx6 = {m: k for k, m in enumerate(rows6)}
    gens = []
    for i in range(6):
        for lm in monomials(1):
            gens.append(pshift(Fi[i], lm, 1, p))
    G = build_mat(gens, rows6, p).transpose()  # rows = generators
    Gr, rk = G.rref()
    assert rk == 36, rk
    piv = []
    for r in range(rk):
        j = 0
        while int(Gr[r, j]) == 0: j += 1
        piv.append(j)
    nonpiv = [j for j in range(len(rows6)) if j not in set(piv)]
    # quadratic form coefficients: for each pair (k<=l) of basis elements, the polynomial w_k w_l H_{ik il} (times 2 if k != l)
    eqs = [dict() for _ in nonpiv]
    for a in range(nvar):
        ia, ma = keep[a]
        for b in range(a, nvar):
            ib, mb = keep[b]
            Pm = pshift(H[ia][ib], tuple(x + y for x, y in zip(ma, mb)), 1 if a == b else 2, p)
            if not Pm: continue
            v = nmod_mat(1, len(rows6), p)
            for e, cc in Pm.items(): v[0, idx6[e]] = cc
            # reduce: v - (v restricted to pivots) * Gr
            vp = nmod_mat(1, rk, [int(v[0, j]) for j in piv], p)
            nf = v - vp * Gr
            for r, j in enumerate(nonpiv):
                cc = int(nf[0, j])
                if cc: eqs[r][(a, b)] = cc
    varnames = [f'c{k}' for k in range(nvar)]
    lines = []
    for eq in eqs:
        if eq: lines.append('+'.join(f'{cc}*c{a}*c{b}' if a != b else f'{cc}*c{a}^2' for (a, b), cc in eq.items()))
    fn = f'/tmp/hodge/results/dense-rational/verify/quad6_{sector}_p{p}.ms'
    open(fn, 'w').write(','.join(varnames) + '\n' + f'{p}\n' + ',\n'.join(lines) + '\n')
    out = subprocess.run(['msolve', '-g', '1', '-t', '2', '-f', fn, '-o', fn + '.g1'], capture_output=True, text=True, timeout=3000)
    txt = open(fn + '.g1').read()
    # parse leading monomials: msolve -g 1 output lists leading monomials of the reduced GB
    import re
    mons = re.findall(r'[c0-9\^\*]+', txt.replace('[', ' ').replace(']', ' ').replace(',', ' '))
    pure = {}
    for m in mons:
        m = m.strip()
        if not m or 'c' not in m: continue
        vs = re.findall(r'c(\d+)(?:\^(\d+))?', m)
        if len(vs) == 1:
            k = int(vs[0][0]); e = int(vs[0][1]) if vs[0][1] else 1
            pure[k] = min(pure.get(k, 99), e)
    allpure = all(k in pure for k in range(nvar))
    results[sector] = {'nvar': nvar, 'neq': len(lines), 'n_leading': len([m for m in mons if 'c' in m]), 'pure_powers': dict(sorted(pure.items())), 'all_vars_have_pure_power': allpure, 'seconds': round(time.time() - t0, 1)}
    print(sector, json.dumps(results[sector]), flush=True)
    with open('/tmp/hodge/results/dense-rational/verify/quad6_v.jsonl', 'a') as f: f.write(json.dumps({'p': p, 'sector': sector, **results[sector]}) + '\n')
ok = all(r['all_vars_have_pure_power'] for r in results.values())
log('excl R-d: degree-6 Hessian relaxation is projectively empty in even/odd/mixed sectors (msolve -g 1 leading ideal m-primary) => deg D >= 7 for all parities',
    'CONFIRMED' if ok else 'REFUTED', f'p={p}: ' + '; '.join(f"{s}: {r['nvar']} vars, {r['neq']} quadrics, pure powers for all vars: {r['all_vars_have_pure_power']}" for s, r in results.items()),
    'independent generator + different prime; logic: homogeneous ideal m-primary mod p => Macaulay matrix full column rank mod p in some degree => full rank over Q => empty over C (rigorous); requires rank_p J_6 = 36 (checked) so the mod-p quadrics are reductions of integral Q-quadrics')
