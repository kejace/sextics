import sys, json; sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
from lib import *
from hodge import symmetry
G = symmetry.diagonal_group(F)
print(G)
n = G.order
chars = {}
for w in G.elements:
    chars[str(w)] = symmetry.residue_character(F, A, q, w, n)
print('order', n, 'generators', G.generators)
# residue characters for generators
log(dict(exp='a3-symmetry', order=n, generators=[list(g) for g in G.generators],
         residue_chars_of_generators=[symmetry.residue_character(F, A, q, g, n) for g in G.generators],
         n_invariant_elements=sum(1 for v in chars.values() if v == 0)))
