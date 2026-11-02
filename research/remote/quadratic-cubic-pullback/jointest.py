import sys; sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core
from flint import nmod_mat
p = core.PRIMES[0]
F = core.reduce_mod(qcp.F, p)
# rank of span of the 21 second partials (a join F=G(l)+H(l') has 9 vanishing mixed partials -> rank <= 12)
second = [core.pderiv(core.pderiv(F,i,p),j,p) for i in range(6) for j in range(i,6)]
mons = core.monomials(4)
col = {e:k for k,e in enumerate(mons)}
flat=[0]*(len(second)*len(mons))
for r,P in enumerate(second):
    for e,v in P.items(): flat[r*len(mons)+col[e]]=v
print('rank of span of 21 second partials:', nmod_mat(len(second),len(mons),flat,p).rank(), '(join would give <= 12)')
# rank of span of the 6 first partials: 6 (nondegenerate)
first=[core.pderiv(F,i,p) for i in range(6)]
mons5=core.monomials(5); col5={e:k for k,e in enumerate(mons5)}
flat=[0]*(6*len(mons5))
for r,P in enumerate(first):
    for e,v in P.items(): flat[r*len(mons5)+col5[e]]=v
print('rank of first partials:', nmod_mat(6,len(mons5),flat,p).rank())
