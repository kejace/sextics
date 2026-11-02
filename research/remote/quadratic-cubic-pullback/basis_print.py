import sys, time
sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
which, Ds = sys.argv[1], sys.argv[2]
if which == 'fermat': F_ = qcp.FERMAT
elif which == 'ours': F_ = qcp.F
else:
    F_, A_, q_, rec = core.load_instance(which)
D = P(Ds)
t0 = time.time()
dim, basis, S_ = fixed_d.target_image(F_, 1, D)
print(f'{which} D={Ds} d={core.degree(D)} imgdim={dim} ({time.time()-t0:.1f}s)')
# print basis elements: leading monomial (pivot) and support size; reduce mod J? just show.
def mstr(e):
    return '*'.join(f'{core.VARS[i]}^{k}' if k>1 else core.VARS[i] for i,k in enumerate(e) if k) or '1'
for b in basis:
    items = sorted(b.items(), reverse=True)
    lead = items[0][0]
    print(f'  lead {mstr(lead):28s} nterms={len(b):3d}  ' + ' '.join(mstr(e) for e,_ in items[1:6]) + (' ...' if len(b)>6 else ''))
