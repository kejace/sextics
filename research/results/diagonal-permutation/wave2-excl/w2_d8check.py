import json, numpy as np, sys
import w2lib as W
from hodge import core, fixed_d
from hodge.core import reduce_mod, sub, add, scale
p = W.P1
D8 = core.parse(json.load(open(f'{W.LANE}/D8_e1kernel.json'))['D'])
red8 = W.Red(8, p, 0, None)
D8p = reduce_mod(D8, p)
nf = red8.nf(D8p)
print('D8 in J?', not np.any(nf), '; NF nonzero coords', int(np.count_nonzero(nf)), 'of', red8.ring_dim, flush=True)
S = W.act_swap01(D8p)
print('(01)D8 - D8 in J?', red8.in_J(sub(S, D8p, p)), '; (01)D8 + D8 in J?', red8.in_J(add(S, D8p, p)), flush=True)
C = W.act_perm(D8p, (1, 2, 0))
M = np.vstack([nf, red8.nf(C)])
print('rank of {NF(D8), NF(cyc D8)}:', W.rank_np(M, p), '(2 = three-cycle moves the line)', flush=True)
print('E1 at P1:', fixed_d.e1_solvable(W.F, W.A, 1, D8, p=p)[0], '; at P2:', fixed_d.e1_solvable(W.F, W.A, 1, D8, p=W.P2)[0], flush=True)
W.log('d8_sector_check', D='D8', d=8, prime=p, result='(01)D8 = +D8 mod J; cyc D8 independent; D8 not in J; E1 holds at P1,P2', note='re-verification of wave-1 d8_probe claims')
