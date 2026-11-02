"""NZD lemma replay at a new prime: for the chart l = v_s + b1 v_1 + b2 v_2 (or l = v_1 + b1 v_2) the quadratic syzygies of the
six restricted partials (7F)_i|_{l=0} vanish for every (b1,b2) in GF(p)-bar  <=>  Ann(coker M(b)^T) = (1)  <=>  (J:l)_6 = J_6 for every such l.
usage: python3 ref_nzd.py"""
import sys, itertools, subprocess, time
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl/referee')
import ref_lib as L
p = L.P1
F, Fi, H = L.Fpolys(None)
def subst_linear(P, s, combo):
    """substitute variable s by - sum_j b_j * var_j; polynomials in 8 variables (6 + b1 + b2): exponent 8-tuples"""
    R = {}
    for e, c in P.items():
        k = e[s]
        base = list(e) + [0, 0]; base[s] = 0
        # expand (-(b1 v1 + b2 v2))^k multinomially
        for ks in itertools.product(range(k + 1), repeat=len(combo)):
            if sum(ks) != k: continue
            coef = c * (-1) ** k
            from math import factorial
            coef = coef * factorial(k)
            for kk in ks: coef //= factorial(kk) if isinstance(coef, int) else 1
            ee = list(base)
            for (bj, vj), kk in zip(combo, ks):
                ee[vj] += kk; ee[6 + bj] += kk
            R[tuple(ee)] = R.get(tuple(ee), 0) + coef
    return {e: c for e, c in R.items() if c}
def chart(name, s, combo, keepvars):
    Fbar = [subst_linear(Fi[i], s, combo) for i in range(6)]
    quads = [m for m in L.monomials(2) if m[s] == 0]           # 15 quadratic monomials in the 5 remaining variables
    cols = [(i, m) for i in range(6) for m in quads]
    rows = {}; ent = {}
    for c, (i, m) in enumerate(cols):
        for e, v in Fbar[i].items():
            x = tuple(e[t] + m[t] for t in range(6)); b = (e[6], e[7])
            r = rows.setdefault(x, len(rows))
            ent[(r, c, b)] = ent.get((r, c, b), 0) + v
    nr, nc = len(rows), len(cols)
    lines = [f'ring r = {p}, (b1,b2), dp;', f'matrix M[{nr}][{nc}];']
    byrc = {}
    for (r, c, b), v in ent.items():
        if v % p: byrc.setdefault((r, c), []).append(f'({v % p})*b1^{b[0]}*b2^{b[1]}')
    for (r, c), terms in byrc.items(): lines.append(f'M[{r+1},{c+1}] = ' + '+'.join(terms) + ';')
    lines += ['module N = module(transpose(M));', 'ideal ann = quotient(N, freemodule(%d));' % nc, 'ideal sann = std(ann);',
              f'"{name}", "rows", {nr}, "cols", {nc}, "ann_size", size(sann), "dim", dim(sann), "ann_is_unit", (size(sann) == 1 and deg(sann[1]) == 0);', 'quit;']
    path = f'nzd_{name}_p{p}.sing'; open(path, 'w').write('\n'.join(lines) + '\n')
    t = time.time(); r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=3000)
    print(r.stdout.strip(), r.stderr.strip()[:200], f'[{time.time()-t:.1f}s]', flush=True)
chart('Xchart0_X0+b1X1+b2X2', 0, [(0, 1), (1, 2)], None)
chart('Xchart1_X1+b1X2', 1, [(0, 2)], None)
chart('Ychart0_Y0+b1Y1+b2Y2', 3, [(0, 4), (1, 5)], None)
chart('Ychart1_Y1+b1Y2', 4, [(0, 5)], None)
print('NZD_DONE (coordinate points X2, Y2 and all six coordinate forms: see referee_colon json, colon_linear)')
