"""Degree-5 denominators: by the E1 kernel computation, any D of degree 5 passing E1 lies in J_5 = span(F_i).
Test structured combinations D = sum c_i F_i (D and D^2), plus the served/new statistic at D."""
import sys, time, itertools, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, served
from hodge import core
from fractions import Fraction
F, A, q, rec = core.load_instance('cycle-q2')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
Fi = [core.pderiv(F, i) for i in range(6)]
t0 = time.time()
def run(c, note):
    D = {}
    for j in range(6):
        if c[j]: D = core.add(D, core.scale(Fi[j], Fraction(c[j])))
    Ds = 'sum c_i F_i, c=' + str(list(c))
    r = served.served_dim(F, A, q, D); served.log_served(D, 'deg5-J', note, r, Dstr=Ds)
    s = f'{note:28s} d=5 |H|={r["order"]} new {r["new"]} (unserved {r["unserved"]}/{r["baseline"]} of {r["sector"]}) cons={r["consistent"]} rk={r["rank"]}'
    if r['consistent']: print('!!!!! CONSISTENT (D) !!!!!', note, flush=True)
    D2 = core.mul(D, D)
    r2 = eq.test_D_fast(F, A, q, D2, 'deg5-J', note=note + ' [D^2]', verbose=False, Dstr='(' + Ds + ')**2')
    rb = eq.test_D_fast(F, A, q, D2, 'deg5-J', note=note + ' [D^2] 2nd prime', verbose=False, Dstr='(' + Ds + ')**2', p=4194287)
    s += f' | D^2 d=10: cons={r2["consistent"]} rk={r2["rank"]} nul={r2["nullity"]} cols={r2["cols"]} |H|={r2["order"]} (p2 cons={rb["consistent"]} rk={rb["rank"]})'
    if r2['consistent'] or rb['consistent']: print('!!!!! CONSISTENT (D^2) !!!!!', note, flush=True)
    print(s + f' [{time.time()-t0:.0f}s]', flush=True)
for i in range(6): run([1 if j == i else 0 for j in range(6)], f'F_{V[i]}')
run([1, 1, 1, 0, 0, 0], 'F_X0+F_X1+F_X2'); run([0, 0, 0, 1, 1, 1], 'F_Y0+F_Y1+F_Y2'); run([1, 1, 1, 1, 1, 1], 'sum all F_i')
for i in range(3): run([1 if j in (i, i + 3) else 0 for j in range(6)], f'F_X{i}+F_Y{i}')
for i in range(3): run([1 if j == i else (-1 if j == i + 3 else 0) for j in range(6)], f'F_X{i}-F_Y{i}')
run([1, -1, 0, 0, 0, 0], 'F_X0-F_X1'); run([1, 0, 0, 0, 1, 0], 'F_X0+F_Y1'); run([1, 0, 0, 0, 0, 1], 'F_X0+F_Y2'); run([0, 1, 0, 0, 0, 1], 'F_X1+F_Y2')
run([1, 2, 3, 0, 0, 0], 'F_X0+2F_X1+3F_X2'); run([1, 1, 1, -1, -1, -1], 'sum F_X - sum F_Y')
random.seed(11)
for k in range(6): run([random.randint(-4, 4) for _ in range(6)], f'random-{k}')
print('DONE', flush=True)
