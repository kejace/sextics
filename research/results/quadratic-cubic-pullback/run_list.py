import sys, time, json
sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
which = sys.argv[1]
F_ = qcp.FERMAT if which == 'fermat' else qcp.F
A_ = qcp.A
items = [l.strip() for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
for line in items:
    fam, Ds = line.split('|', 1)
    fam = fam.strip(); Ds = Ds.strip()
    D = eval(Ds, {'P': P, 'Q': qcp.Q, 'core': core, 'delta': P(open('delta.txt').read()), 'FERMAT': qcp.FERMAT, 'F': qcp.F})
    if isinstance(D, str): D = P(D)
    t0 = time.time()
    dim, basis, S_, dt = qcp.image_dim(D, F_=F_)
    r = qcp.test(D, fam + ('/fermat' if which=='fermat' else ''), note=f'{Ds} | image dim {dim}', F_=F_, A_=A_)
    print(f"{which:6s} {fam:18s} d={core.degree(D)} {Ds[:60]:60s} E1={r['e1']} cons={r.get('consistent')} rank={r.get('rank')} null={r.get('nullity')} imgdim={dim} {time.time()-t0:.1f}s", flush=True)
