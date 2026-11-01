import sys, time, json
sys.path.insert(0,'/tmp/hodge/results/quadratic-cubic-pullback')
import qcp
from qcp import core, fixed_d, P
which = sys.argv[1]
F_ = qcp.FERMAT if which == 'fermat' else qcp.F
A_ = qcp.A
items = [l.strip() for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
delta = P(open('/tmp/hodge/results/quadratic-cubic-pullback/delta.txt').read())
def Fpart(xdeg):
    return {e:c for e,c in F_.items() if e[0]+e[1]+e[2]==xdeg}
def EI(I, G=None):
    G = F_ if G is None else G
    out = {}
    for i in I:
        out = core.add(out, core.mul(qcp.var(i), core.pderiv(G, i)))
    return out
env = {'P': P, 'Q': qcp.Q, 'core': core, 'delta': delta, 'FERMAT': qcp.FERMAT, 'F': F_, 'Fpart': Fpart, 'EI': EI, 'Fraction': qcp.Fraction}
for line in items:
    parts = line.split('|')
    fam, Ds = parts[0].strip(), parts[1].strip()
    flags = parts[2].strip() if len(parts) > 2 else ''
    D = eval(Ds, env)
    if isinstance(D, str): D = P(D)
    t0 = time.time()
    dim = None
    if 'img' in flags:
        dim, basis, S_, dt = qcp.image_dim(D, F_=F_)
    r = qcp.test(D, fam + ('/fermat' if which=='fermat' else ''), note=f'{Ds} | image dim {dim}', F_=F_, A_=A_)
    print(f"{which:6s} {fam:14s} d={core.degree(D)} {Ds[:70]:70s} E1={r['e1']} cons={r.get('consistent')} rank={r.get('rank')} null={r.get('nullity')} imgdim={dim} {time.time()-t0:.1f}s", flush=True)
    if r.get('consistent') and r.get('N') is not None:
        qcp.save_witness(D, r['N'], r['Q'], path=f'/tmp/hodge/results/quadratic-cubic-pullback/witness-candidate-{which}-{abs(hash(Ds))%100000}.json')
        print('   >>> WITNESS SAVED', flush=True)
