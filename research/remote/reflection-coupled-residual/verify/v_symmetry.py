"""Independent check of the diagonal/monomial symmetry group of F and residue characters."""
import sys, itertools, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
Fs, As, qv = load_rec(); F = parse(Fs); A = parse(As)
supp = sorted(F)
# Any diagonal projective symmetry diag(z_i) with F(z Z) = lam F: Z_i^6 in F forces z_i^6 = lam for all i;
# normalize z_0 = 1 -> lam = 1, z_i in mu_6.  Exhaustive search over weights w_i in Z/6 (w_0 = 0).
sols = []
for w in itertools.product(range(6), repeat=5):
    w = (0,) + w
    if all(sum(wi*ei for wi, ei in zip(w, e)) % 6 == 0 for e in supp):
        sols.append(w)
print('diagonal projective symmetries (weights mod 6, w0=0):', sols)
# permutation symmetries (support-preserving permutations; a monomial symmetry perm*diag needs the permuted support = support)
perms = []
for pi in itertools.permutations(range(6)):
    if sorted(tuple(e[pi[i]] for i in range(6)) for e in supp) == supp:
        perms.append(pi)
print('support-preserving permutations:', perms)
# X<->Y exchange
sw = sorted(e[3:]+e[:3] for e in supp)
print('F(Y,X)==F as support?', sw == supp)
# residue characters of A = X0^3 Y0^3 and companion X0^3 Y2^3 under r, rho, tau (sign actions): numerator sign * det
acts = {'r': (1,1,1,1,1,-1), 'rho': (1,1,1,-1,-1,1), 'tau': (1,1,1,-1,-1,-1)}
def numsign(e, s):
    v = 1
    for ei, si in zip(e, s):
        if si == -1 and ei % 2: v = -v
    return v
for name, s in acts.items():
    det = 1
    for si in s: det *= si
    Fsign = {numsign(e, s) for e in supp}
    for Aexpr in ('X0**3*Y0**3', 'X0**3*Y2**3'):
        e = next(iter(parse(Aexpr)))
        print(f'{name}: F sign {Fsign}, det {det}, A={Aexpr} numerator sign {numsign(e,s)}, residue char {numsign(e,s)*det}')
# parity facts about F
print('F even in each variable individually?', [all(e[i] % 2 == 0 for e in supp) for i in range(6)])
print('F Y-even?', all((e[3]+e[4]+e[5]) % 2 == 0 for e in supp), ' F Y2-even?', all(e[5] % 2 == 0 for e in supp))
# smoothness-at-5 argument: 5/7 mod 5
print('5/7 mod 5 =', (5 * pow(7, -1, 5)) % 5)
logjson(VLOG, {'check': 'symmetry', 'diag_weights_mod6': sols, 'perms': [list(p) for p in perms], 'XY_swap': sw == supp})
