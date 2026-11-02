"""Independent rebuild of the degree-6 Hessian relaxation: D = w(F), w linear vector field (36 coeffs; even part
18, odd part 18).  Necessary condition (derived in the verify note): Phi(w) = w^T Hess(F) w has even part in J_6^+
and odd part in J_6^-.  Build the quadratic equations mod p (and over Q for pure parities) and solve with Singular.
Usage: relax6.py <p|0> <e|o|m>"""
import sys, json, subprocess, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation/verify')
from vlib import *
p = int(sys.argv[1]); mode = sys.argv[2]
Fs, As, q = sympy_F()
F = poly_from_sympy(Fs)
Fq = F if p == 0 else reduce_mod(F, p)
Fi = [pderiv(Fq, i, p or None) for i in range(6)]
Hess = [[pderiv(Fi[i], j, p or None) for j in range(6)] for i in range(6)]
# unknown vector fields: component j gets linear monomials of parity (parD - par(Z_j))
def basis(parD):
    out = []
    for j in range(6):
        pj = 0 if j < 3 else 1
        for m in monomials(1, (parD - pj) % 2):
            out.append((j, m))
    return out
B = basis(0) if mode == 'e' else basis(1) if mode == 'o' else basis(0) + basis(1)
nb = len(B)
# Phi(w) = sum_{a,b} w_a w_b m_a m_b Hess[j_a][j_b]; collect per quadratic monomial a<=b
terms = {}
for a in range(nb):
    for b in range(a, nb):
        (j, m), (l, m2) = B[a], B[b]
        mm = tuple(x + y for x, y in zip(m, m2))
        P = shift(Hess[j][l], mm, 1 if a == b else 2, p or None)
        if P: terms[(a, b)] = P
# split by parity and reduce modulo J_6^{par}
import sympy
eqs = []
for par in (0, 1):
    mons = monomials(6, par); idx = {m: i for i, m in enumerate(mons)}
    gens = jac_gens(Fq, 6, par, p or None) if p else None
    if p:
        from flint import nmod_mat
        K = len(gens); n = len(mons)
        flat = [0] * (K * n)
        for k, g in enumerate(gens):
            for m, v in g.items(): flat[k * n + idx[m]] = v
        Mj = nmod_mat(K, n, flat, p)
        Rr0, r = Mj.rref()
        piv = []
        for i in range(r):
            j = 0
            while int(Rr0[i, j]) == 0: j += 1
            piv.append(j)
        Rr = nmod_mat(r, n, [int(Rr0[i, j]) for i in range(r) for j in range(n)], p)
        keys = [k for k, P in terms.items() if parity_of(P) == par]
        if not keys: continue
        flat = []
        for k in keys:
            row = [0] * n
            for m, v in terms[k].items(): row[idx[m]] = v
            flat.extend(row)
        U = nmod_mat(len(keys), n, flat, p)
        Up = nmod_mat(len(keys), r, [int(U[i, pj]) for i in range(len(keys)) for pj in piv], p)
        Ured = U - Up * Rr
        nonpiv = [j for j in range(n) if j not in set(piv)]
        for col in nonpiv:
            eq = {k: int(Ured[i, col]) for i, k in enumerate(keys) if int(Ured[i, col])}
            if eq: eqs.append(eq)
        print(f'par {par}: dim S_6 = {n}, rank J_6 = {r} (expected {hilb()[6]["J_plus" if par == 0 else "J_minus"]}), eqs so far {len(eqs)}', flush=True)
    else:
        # exact over Q: use sympy Matrix rref on the J generators (18 x n) -- small
        gens = jac_gens(F, 6, par, None)
        K = len(gens); n = len(mons)
        Mj = sympy.Matrix(K, n, lambda k, i: gens[k].get(mons[i], 0))
        Rr, piv = Mj.rref()
        r = len(piv)
        keys = [k for k, P in terms.items() if parity_of(P) == par]
        if not keys: continue
        U = sympy.Matrix(len(keys), n, lambda i, j: terms[keys[i]].get(mons[j], 0))
        Up = U[:, list(piv)]
        Ured = U - Up * Rr[:r, :]
        nonpiv = [j for j in range(n) if j not in set(piv)]
        for col in nonpiv:
            eq = {k: Ured[i, col] for i, k in enumerate(keys) if Ured[i, col] != 0}
            if eq: eqs.append(eq)
        print(f'par {par}: dim S_6 = {n}, rank J_6 = {r}, eqs so far {len(eqs)}', flush=True)
path = f'relax6_{mode}_p{p}.sing'
with open(path, 'w') as f:
    f.write(f'ring R = {p}, (a(1..{nb})), dp;\nideal I;\n')
    for eq in eqs:
        f.write('I[size(I)+1] = ' + '+'.join(f'({v})*a({a+1})*a({b+1})' for (a, b), v in eq.items()) + ';\n')
    f.write('ideal G = std(I);\nprint("dim = " + string(dim(G)));\nprint("mult = " + string(mult(G)));\n')
    f.write('LIB "primdec.lib"; list Lp = minAssGTZ(G); print("ncomp = " + string(size(Lp))); Lp;\n')
    f.write('ideal Rad = radical(G); print("radical:"); Rad;\nquit;\n')
print('unknowns', nb, 'equations', len(eqs), '->', path, flush=True)
t0 = time.time()
out = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=7200).stdout
print(out, flush=True)
print(f'Singular time {time.time()-t0:.0f}s')
logjson(claim=f'degree-6 Hessian relaxation, mode={mode}, p={p}: solution set of Phi(w) in J_6', verdict='DATA', unknowns=nb, equations=len(eqs), singular=out[:2500])
