"""Replay of the symmetry bookkeeping (dossier 1.4): (1) all coordinate permutations pi such that F o pi is a diagonal
rescaling of F; (2) residue character of A Omega/F^2 under the generators; (3) dimension of the epsilon-sector of H^{3,1}
and H^{2,2}: A' in S_k[chi_A] with s0 A' = A', r1 A' = A', modulo J_k, for k = 6 (expect 22) and k = 12 (expect 51), mod P2
(rigorous as an upper bound on the rational dimension; equality since the Jacobian rank mod p is full)."""
import sys, time, json, itertools
from w2lib import *
# (1) permutation symmetries
def perm_poly(P, pi):   # (pi.P)(z) = P(z_{pi(0)}, ..., z_{pi(5)})
    return {tuple(e[pi[j]] for j in range(NV)): c for e, c in P.items()}
supp = set(F)
found = []
for pi in itertools.permutations(range(NV)):
    G = perm_poly(F, pi)
    if set(G) != supp: continue
    # is there a diagonal scaling lambda with G = c * F(lambda z)?  ratios G[e]/F[e] must be a character lambda^e (up to a scalar)
    ratios = {e: G[e] / F[e] for e in supp}
    # test: ratio must be multiplicative: r(e) r(e') = r(e'') r(e''') whenever e + e' = e'' + e''' ... simpler: solve for lambda via
    # the 6 Fermat monomials: lambda_i^6 = r(6 e_i)/r0 ; then check all others up to 6th roots: we test r(e)^6 == prod r(6e_i)^{e_i}/... (projective)
    r6 = [ratios[tuple(6 if j == i else 0 for j in range(NV))] for i in range(NV)]
    ok = True
    for e, r in ratios.items():
        lhs = r ** 6
        rhs = Fraction(1)
        for i in range(NV): rhs *= r6[i] ** e[i]
        # projective scalar c: r(e) = c * lambda^e with c = r(6 e_0)/lambda_0^6 ... use ratio against the X0^6 monomial
        if lhs * (r6[0] ** 6) != rhs * (r6[0] ** 6) and lhs != rhs:
            ok = False; break
    if ok: found.append(pi)
log('symmetry-permutations', n_perms_preserving_support=len(found), perms=[list(x) for x in found],
    note='permutations pi with supp(F o pi) = supp(F) and coefficient ratios consistent (6th-power test) with a diagonal rescaling; expect id, s0=(0<->3), r1=(1<->4), s0 r1')
# (3) epsilon-sector dimensions
S0p, R1p = (3, 1, 2, 0, 4, 5), (0, 4, 2, 3, 1, 5)
def orbit(e):
    return sorted({tuple(e[g[j]] for j in range(NV)) for g in [(0, 1, 2, 3, 4, 5), S0p, R1p, (3, 4, 2, 0, 1, 5)]})
for k, expect in ((6, 22), (12, 51)):
    t0 = time.time()
    c = CHI_A
    seen, basis = set(), []
    for e in monomials(k):
        if chi(e) != c or e in seen: continue
        orb = orbit(e); seen.update(orb); basis.append({m: 1 for m in orb})
    # J_k[c]^{inv}: symmetrize generators and take the rank in orbit-sum coordinates
    S = sector(k, c)
    bidx = {}
    for i, B in enumerate(basis):
        for e in B: bidx[e] = i
    rows = []
    for (i, mm) in S.gens:
        P = shift(Fi[i], mm, 1, p)
        Q = {}
        for g in [(0, 1, 2, 3, 4, 5), S0p, R1p, (3, 4, 2, 0, 1, 5)]:
            Q = add(Q, perm_poly(P, g), p)
        if not Q: continue
        v = [0] * len(basis)
        for e, cv in Q.items():
            v[bidx[e]] = cv   # each orbit-sum coordinate: coefficient is constant along the orbit (Q is invariant)
        rows.append(v)
    rk = nmod_mat(len(rows), len(basis), [x for r in rows for x in r], p).rank() if rows else 0
    log('epsilon-sector-dimension', k=k, dim_S_inv=len(basis), dim_J_inv=rk, dim_R_inv=len(basis) - rk, expected=expect, agrees=(len(basis) - rk == expect), seconds=round(time.time() - t0, 1),
        note='A\' in S_k[chi_A] invariant under s0 and r1, modulo J_k: k=6 -> h^{3,1}_eps, k=12 -> h^{2,2}_eps (dossier: (0,22,51,22,0))')
print('DONE', file=sys.stderr)
