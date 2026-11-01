"""Experiment 3: denominators in the Jacobian ideal (pass E1 trivially): F_i (d=5), x_j F_i (d=6),
F_i F_j (d=10), radial combinations sum a_j x_j F_j (d=6), and squares of the degree-5/6 ones."""
import sys, time, json, itertools, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
from fractions import Fraction

F, A, q, rec = core.load_instance('cycle-q2')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
Fi = [core.pderiv(F, i) for i in range(6)]
x = [core.parse(v) for v in V]
t0 = time.time()
def run(D, fam, note, square=True):
    r = eq.test_D(F, A, q, D, fam, note=note + ' [D]', verbose=False, do_e1=True)
    s = f'{fam:12s} {note[:60]:60s} | d={r["d"]}: e1={r.get("e1")} cons={r.get("consistent")} rk={r.get("rank")} nul={r.get("nullity")} cols={r.get("cols")} |H|={r.get("order")}'
    if square:
        D2 = core.mul(D, D)
        r2 = eq.test_D(F, A, q, D2, fam, note=note + ' [D^2]', verbose=False, do_e1=False, Dstr='(' + core.to_str(D) + ')**2')
        s += f' | D^2 d={r2["d"]}: cons={r2.get("consistent")} rk={r2.get("rank")} nul={r2.get("nullity")} cols={r2.get("cols")}'
        if r2.get('consistent'): print('!!!!! CONSISTENT !!!!!', note, flush=True)
    if r.get('consistent'): print('!!!!! CONSISTENT !!!!!', note, flush=True)
    print(s + f' [{time.time()-t0:.0f}s]', flush=True)

for i in range(6):
    run(Fi[i], 'partial', f'F_{V[i]}')
for i in range(6):
    for j in range(6):
        run(core.mul(x[j], Fi[i]), 'x_j*F_i', f'{V[j]}*F_{V[i]}')
# radial-type combinations (semi-invariant since all x_j F_j are invariant)
combos = {
 'B6=sumX_jF_Xj (excluded by fqlx, control)': [1, 1, 1, 0, 0, 0],
 'sumY_jF_Yj': [0, 0, 0, 1, 1, 1],
 'Grad0 (excluded by fqlx, control)': [1, 0, 0, 1, 0, 0],
 'Grad1': [0, 1, 0, 0, 1, 0],
 'Grad2': [0, 0, 1, 0, 0, 1],
 'X0F0-X1F1': [1, -1, 0, 0, 0, 0],
 'X0F0-Y0F3': [1, 0, 0, -1, 0, 0],
 'X0F0+X1F1-2X2F2': [1, 1, -2, 0, 0, 0],
 'X0F0+Y1F4': [1, 0, 0, 0, 1, 0],
 'X0F0+Y2F5': [1, 0, 0, 0, 0, 1],
 'X0F0+X1F1+Y0F3+Y1F4': [1, 1, 0, 1, 1, 0],
 'X0F0+2X1F1+3X2F2': [1, 2, 3, 0, 0, 0],
 'X0F0-X1F1+X2F2-Y0F3+Y1F4-Y2F5': [1, -1, 1, -1, 1, -1],
}
random.seed(7)
for k in range(4):
    combos[f'random-{k}'] = [random.randint(-5, 5) for _ in range(6)]
for name, a in combos.items():
    D = {}
    for j in range(6):
        if a[j]:
            D = core.add(D, core.scale(core.mul(x[j], Fi[j]), Fraction(a[j])))
    if not D: continue
    run(D, 'radial', name)
# products of partials (d=10) -- no square (d=20 too big here)
for i in range(6):
    for j in range(i, 6):
        run(core.mul(Fi[i], Fi[j]), 'F_i*F_j', f'F_{V[i]}*F_{V[j]}', square=False)
print('DONE', flush=True)
