import sys; sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
from rcr import *
import time
t0=time.time()
J = get_J6()
print('S_6^- dim', len(J.mons), 'J_6^- dim (rank)', J.rank, 'R_6^- dim', J.quot_dim, f'{time.time()-t0:.1f}s')
ev = J.sector_mask(0); od = J.sector_mask(1)
print('R_6^- split by Y2-parity of nonpivot monomials: even-Y2 (anti-invariant residue sector):', len(ev), ' odd-Y2 (invariant, descends to W):', len(od))
# cross-check: Fermat count
import itertools
cnt_e=cnt_o=0
for m in monomials(6,1):
    if max(m)<=4:
        if m[5]%2==0: cnt_e+=1
        else: cnt_o+=1
print('Fermat monomial-basis count even/odd Y2:', cnt_e, cnt_o)
av = A_vec()
print('A normal form nonzero coords:', sum(1 for x in av if x), 'in even sector?', all(av[k]==0 for k in od))
