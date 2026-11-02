import sys, time, json, itertools
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
from fractions import Fraction
F, A, q, rec = core.load_instance('cyclic-triple-residual')
print('q =', q, ' #F terms', len(F), ' A =', core.to_str(A))
# diagonal symmetry group: gamma_i = zeta6^{a_i}
supp = list(F.keys())
sols = []
for a in itertools.product(range(6), repeat=6):
    ok = True
    e0 = supp[0]
    w0 = sum(ai*ei for ai, ei in zip(a, e0)) % 6
    for e in supp:
        if sum(ai*ei for ai, ei in zip(a, e)) % 6 != w0:
            ok = False; break
    if ok: sols.append((a, w0))
print('|Gamma| (incl. scalars) =', len(sols), ' projective order =', len(sols)//6)
# residue character of A Omega/F^2: chi_A * det / mu
Ae = next(iter(A.keys()))
chars = {}
for a, w0 in sols:
    chiA = sum(ai*ei for ai, ei in zip(a, Ae)) % 6
    det = sum(a) % 6
    c = (chiA + det - w0) % 6
    chars[a] = c
# projective classes: mod scalars
seen = set(); reps = []
for a, w0 in sols:
    key = tuple((ai - a[0]) % 6 for ai in a)
    if key in seen: continue
    seen.add(key); reps.append((key, chars[tuple(key)] if tuple(key) in chars else None))
print('projective elements (normalized a0=0) and residue character exponent (zeta6^c):')
for key, c in reps:
    print('  ', key, ' char exp', chars[key])
inv = [k for k, c in reps if chars[k] == 0]
print('elements fixing the residue class:', inv)
# junior elements (age < 1) among fixed points: compute age on P^5 at fixed point components.
# For diagonal gamma, fixed loci are coordinate subspaces where the eigenvalue is constant.
def age_info(key):
    # eigenvalue exponents b_i = a_i mod 6, projective: shift all by s to normalize each fixed component
    out = []
    for s in set(key):
        comp = [i for i in range(6) if key[i] == s]
        # fixed component P(comp); normal eigenvalues (a_i - s) mod 6 for i not in comp
        norm = [((key[i]-s) % 6) for i in range(6) if i not in comp]
        age = sum(norm)/6.0
        # restricted to X: tangent space of X at a fixed point drops one eigen-direction: the normal of X in P^5 has eigenvalue mu*(scalar)...
        out.append((comp, norm, age))
    return out
print('ages on P^5 fixed components (comp, normal eigen-exponents, age=sum/6):')
for key, c in reps:
    if key == (0,)*6: continue
    print('  ', key, 'char', chars[key], age_info(key))
