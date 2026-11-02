"""Re-verify the representation-theoretic claims of the dossier (§1.4-1.5) with independent code:
(i) projective diagonal symmetry group of F (engine symmetry.diagonal_group), (ii) the signed-permutation symmetries of F,
(iii) S3-character of R_6^- (primitive (3,1), dim 208) => multiplicities (triv, sign, std), and the (01)-eigenspace dims,
(iv) the same for R_12^- and R_6^+.  Logs exp='reptheory'."""
import sys, time, itertools, json
import numpy as np
import w2lib as W
from hodge import symmetry, core
from hodge.core import reduce_mod, monomials
p = W.P1
F = W.F
t0 = time.time()
G = symmetry.diagonal_group(F)
res = {'diagonal_group_order': G.order if hasattr(G, 'order') else None}
# signed permutations preserving F up to scalar
perms = list(itertools.permutations(range(6)))
keep = []
Fset = F
for perm in perms:
    for signs in itertools.product((1, -1), repeat=6):
        P = {}
        for e, c in F.items():
            ne = [0] * 6; sg = 1
            for i in range(6):
                ne[perm[i]] += e[i]
                if e[i] % 2 == 1 and signs[i] == -1: sg = -sg
            P[tuple(ne)] = P.get(tuple(ne), 0) + sg * c
        P = {e: c for e, c in P.items() if c}
        # proportional to F?
        if set(P) == set(F):
            ratios = {P[e] / F[e] for e in F}
            if len(ratios) == 1:
                keep.append((perm, signs))
res['signed_permutation_symmetries'] = len(keep)
res['signed_permutation_list'] = [(pm, sg) for pm, sg in keep][:24]
# characters on R_n^par
def char_table(n, par):
    red = W.Red(n, p, par, None)
    std = red.std_mons
    sidx = {m: i for i, m in enumerate(std)}
    def trace_of(perm):
        # action: X_i -> X_perm[i], Y_i -> Y_perm[i] on monomials; NF; trace of the matrix on standard monomials
        imgs = []
        for m in std:
            ne = [0] * 6
            for i in range(3):
                ne[perm[i]] += m[i]; ne[3 + perm[i]] += m[3 + i]
            imgs.append({tuple(ne): 1})
        M = red.nf_polys(imgs)
        return int(np.trace(M) % p)
    chi_e = len(std)
    chi_t = trace_of((1, 0, 2))
    chi_c = trace_of((1, 2, 0))
    # symmetric values (traces are small integers)
    sym = lambda v: v if v <= p // 2 else v - p
    chi_t, chi_c = sym(chi_t), sym(chi_c)
    triv = (chi_e + 3 * chi_t + 2 * chi_c) // 6
    sign = (chi_e - 3 * chi_t + 2 * chi_c) // 6
    std_mult = (2 * chi_e - 2 * chi_c) // 6
    return {'dim': chi_e, 'chi_transposition': chi_t, 'chi_3cycle': chi_c, 'triv': triv, 'sign': sign, 'std_mult': std_mult,
            'check': triv + sign + 2 * std_mult == chi_e, 'dim_01_plus': (chi_e + chi_t) // 2, 'dim_01_minus': (chi_e - chi_t) // 2}
for (n, par, label) in [(6, 1, 'R6_odd (primitive (3,1) tau-invariant)'), (12, 1, 'R12_odd'), (6, 0, 'R6_even'), (8, 0, 'R8_even'), (8, 1, 'R8_odd'), (9, 0, 'R9_even'), (10, 0, 'R10_even'), (10, 1, 'R10_odd'), (16, 1, 'R16_odd'), (17, 1, 'R17_odd')]:
    try:
        res[label] = char_table(n, par)
        print(label, res[label], flush=True)
    except Exception as ex:
        res[label] = {'error': str(ex)}
res['seconds'] = round(time.time() - t0, 1)
json.dump(res, open(f'{W.PKT}/reptheory.json', 'w'), indent=1, default=str)
W.log('reptheory', prime=p, result=f"diag group order {res['diagonal_group_order']}; signed perms {len(keep)}; R6^-: {res['R6_odd (primitive (3,1) tau-invariant)']}", seconds=res['seconds'])
print(json.dumps({k: v for k, v in res.items() if k != 'signed_permutation_list'}, default=str)[:1500])
