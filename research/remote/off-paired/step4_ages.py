"""Age criterion: enumerate the diagonal symmetry group G of F (mod scalars), its eigenspaces on P^5,
fixed components on X, and tangent ages (Kollar-Larsen / book sec 3.6)."""
import itertools
from fractions import Fraction
V = ('X0','X1','X2','Y0','Y1','Y2')
terms = [((3,1,0,1,1,0)), ((0,3,1,0,1,1)), ((2,1,1,1,0,1))]   # exponents of the mixed monomials of F
G = []
for w in itertools.product(range(6), repeat=6):
    if w[0] != 0: continue   # normalize modulo scalars
    if all(sum(a*b for a, b in zip(t, w)) % 6 == 0 for t in terms):
        G.append(w)
print('|G/scalars| =', len(G))
# residue character check: A=X0^3Y0^3, Omega weight = sum w
bad = [w for w in G if (3*w[0] + 3*w[3] + sum(w)) % 6 != 0]
print('elements under which the residue is NOT invariant:', bad)
rows = []
for w in G:
    if all(x == 0 for x in w): continue
    eig = {}
    for i, x in enumerate(w): eig.setdefault(x, []).append(i)
    comps = []
    for lam, idx in eig.items():
        if len(idx) < 2: continue
        age = Fraction(sum(((w[j] - lam) % 6) for j in range(6) if j not in idx), 6)
        comps.append((lam, [V[j] for j in idx], age))
    minage = min((c[2] for c in comps), default=None)
    rows.append((w, comps, minage))
rows.sort(key=lambda r: (r[2] is None, r[2]))
for w, comps, minage in rows:
    print(w, 'min age', minage, ' components:', [(lam, ''.join(v), str(a)) for lam, v, a in comps])
junior = [r for r in rows if r[2] is not None and 0 < r[2] < 1]
print('\nJUNIOR (0<age<1) elements:', [(r[0], str(r[2])) for r in junior])
