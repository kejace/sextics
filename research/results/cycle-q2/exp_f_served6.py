"""Served-dimension survey at D itself (degree <= 6, cheap): which natural denominators reach any of the
90-dim (A-character) degree-12 target sector."""
import sys, time, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, served
from hodge import core
from fractions import Fraction
F, A, q, rec = core.load_instance('cycle-q2')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
t0 = time.time()
def run(Dstr, fam, note=''):
    D = core.parse(Dstr)
    r = served.served_dim(F, A, q, D)
    served.log_served(D, fam, note or Dstr, r)
    print(f'{fam:16s} {Dstr[:80]:80s} d={core.degree(D)} |H|={r["order"]} served {r["served"]}/{r["sector"]} cons={r["consistent"]} rk={r["rank_NQ"]}/{r["cols_NQ"]} [{time.time()-t0:.0f}s]', flush=True)
    return r
run('1', 'sanity')
run('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6', 'pure-power')
run('X0**3*X1*Y0*Y1 + X1**3*X2*Y1*Y2 + X0*X2**3*Y0*Y2', 'coupling-sum')
run('X0**6 + X1**6 + X2**6', 'X-block'); run('Y0**6 + Y1**6 + Y2**6', 'Y-block')
for i in range(3): run(f'X{i}**6 + Y{i}**6', 'pair-fermat')
run('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2', 'edge'); run('X1**6 + X2**6 + Y1**6 + Y2**6 + X1**3*X2*Y1*Y2/2', 'edge'); run('X0**6 + X2**6 + Y0**6 + Y2**6 + X0*X2**3*Y0*Y2/2', 'edge')
for c in ['X0**3*X1*Y0*Y1', 'X1**3*X2*Y1*Y2', 'X0*X2**3*Y0*Y2']: run(c, 'single-coupling')
for quad in itertools.combinations(range(6), 4): run(' + '.join(f'{V[i]}**6' for i in quad), 'fermat-cone4')
for tri in itertools.combinations(range(6), 3): run(' + '.join(f'{V[i]}**3' for i in tri), 'fermat-cone3')
for tri in itertools.combinations(range(6), 3): run(' + '.join(f'{V[i]}**6' for i in tri), 'fermat-cone3-d6')
for tri in itertools.combinations(range(6), 3): run(' + '.join(f'{V[i]}**4' for i in tri), 'fermat-cone3-d4')
Fi = [core.pderiv(F, i) for i in range(6)]
for i in range(6):
    r = served.served_dim(F, A, q, Fi[i]); served.log_served(Fi[i], 'partial', f'F_{V[i]}', r)
    print(f'partial F_{V[i]}: served {r["served"]}/{r["sector"]} cons={r["consistent"]}', flush=True)
for i in range(6):
    D = core.mul(core.parse(V[i]), Fi[i]); r = served.served_dim(F, A, q, D); served.log_served(D, 'radial', f'{V[i]}F_{V[i]}', r)
    print(f'radial {V[i]}F_{V[i]}: served {r["served"]}/{r["sector"]} cons={r["consistent"]}', flush=True)
for a, b, c in [('X0**4', '0', 'Y0**4'), ('X0**4 + X2**4', '0', 'Y0**4 + Y2**4'), ('X2**4', '0', 'Y2**4'), ('X0**4', '0', 'X2**4'),
                ('Y0**4', '0', 'Y2**4'), ('X0**2*X2**2', 'X0*X2*Y0*Y2', 'Y0**2*Y2**2'), ('X0**4 + Y0**4', '0', 'X2**4 + Y2**4'),
                ('X0**4', 'X0*X2*Y0*Y2', 'Y0**4'), ('X0**3*Y0', '0', 'X2**3*Y2'), ('X0*X2**3', '0', 'Y0*Y2**3')]:
    Dstr = f'X1**2*({a}) + X1*Y1*({b}) + Y1**2*({c})'.replace(' + X1*Y1*(0)', '')
    run(Dstr, 'double-P3')
print('DONE', flush=True)
