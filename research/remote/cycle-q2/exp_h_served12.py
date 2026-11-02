"""D^2 served/new statistic (degree 12) for the most natural sextic supports."""
import sys, time, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, served
from hodge import core
F, A, q, rec = core.load_instance('cycle-q2')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
t0 = time.time()
def run(Dstr, fam, note=''):
    D = core.parse(Dstr); D2 = core.mul(D, D)
    r = served.served_dim(F, A, q, D2, with_consistency=False)
    served.log_served(D2, fam, (note or Dstr) + ' [D^2]', r, Dstr='(' + Dstr + ')**2')
    print(f'{fam:16s} ({Dstr[:70]})^2 |H|={r["order"]} new {r["new"]} (unserved {r["unserved"]}/{r["baseline"]} of {r["sector"]}) rkNQ={r["rank_NQ"]}/{r["cols_NQ"]} [{time.time()-t0:.0f}s]', flush=True)
which = sys.argv[1] if len(sys.argv) > 1 else 'all'
if which in ('all', 'a'):
    run('X0**6 + X1**6 + X2**6 + Y1**6', 'fermat-cone4')
    run('X0**6 + X1**6 + Y0**6 + Y1**6', 'fermat-cone4')
    run('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2', 'edge-cone')
    run('X1**6 + X2**6 + Y1**6 + Y2**6 + X1**3*X2*Y1*Y2/2', 'edge-cone')
    run('X0**6 + X2**6 + Y0**6 + Y2**6 + X0*X2**3*Y0*Y2/2', 'edge-cone')
    run('X0**6 + X1**6 + X2**6', 'X-block')
    run('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6', 'pure-power')
    run('X0**3*X1*Y0*Y1 + X1**3*X2*Y1*Y2 + X0*X2**3*Y0*Y2', 'coupling-sum')
if which in ('all', 'b'):
    for quad in itertools.combinations(range(6), 4):
        if quad in ((0, 1, 2, 4), (0, 1, 3, 4)): continue
        run(' + '.join(f'{V[i]}**6' for i in quad), 'fermat-cone4')
print('DONE', flush=True)
