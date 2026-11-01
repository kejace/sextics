"""Worker B: Fermat cones (4-variable sextic cones with line vertex), named F-subset supports (D and D^2),
Jacobian-ideal family (F_i, x_j F_i, radial combos, F_i F_j).  Fast modular path."""
import sys, time, itertools, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq
from hodge import core
from fractions import Fraction
F, A, q, rec = core.load_instance('cycle-q2')
V = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
t0 = time.time()
def run(D, fam, note, square=True, e1=True):
    Ds = core.to_str(D)
    r = eq.test_D_fast(F, A, q, D, fam, note=note + ' [D]', verbose=False, do_e1=e1)
    s = f'{fam:14s} {note[:70]:70s} | d={r["d"]}: e1={r.get("e1")} cons={r.get("consistent")} rk={r.get("rank")} nul={r.get("nullity")} cols={r.get("cols")} |H|={r.get("order")}'
    if r.get('consistent'): print('!!!!! CONSISTENT (D) !!!!!', note, Ds, flush=True)
    if square:
        D2 = core.mul(D, D)
        r2 = eq.test_D_fast(F, A, q, D2, fam, note=note + ' [D^2]', verbose=False, do_e1=False, Dstr='(' + Ds + ')**2')
        s += f' | D^2 d={r2["d"]}: cons={r2.get("consistent")} rk={r2.get("rank")} nul={r2.get("nullity")} cols={r2.get("cols")} |H|={r2.get("order")}'
        if r2.get('consistent'): print('!!!!! CONSISTENT (D^2) !!!!!', note, Ds, flush=True)
    print(s + f' [{time.time()-t0:.0f}s]', flush=True)

which = sys.argv[1] if len(sys.argv) > 1 else 'all'
if which in ('all', 'cones'):
    # sextic Fermat cones on 4 variables (line vertex), with +/- sign variants on one term
    for quad in itertools.combinations(range(6), 4):
        D = core.parse(' + '.join(f'{V[i]}**6' for i in quad))
        run(D, 'fermat-cone4', 'cone ' + ','.join(V[i] for i in quad))
    # cubic cones in 3 variables (plane vertex): Fermat cubics x^3+y^3+z^3 on all triples, D^2 degree 6 (cheap)
    for tri in itertools.combinations(range(6), 3):
        D = core.parse(' + '.join(f'{V[i]}**3' for i in tri))
        run(D, 'fermat-cone3', 'cubic cone ' + ','.join(V[i] for i in tri))
    # quartic/quintic/sextic Fermat cones in 3 variables
    for dd in [4, 5, 6]:
        for tri in itertools.combinations(range(6), 3):
            D = core.parse(' + '.join(f'{V[i]}**{dd}' for i in tri))
            run(D, f'fermat-cone3-d{dd}', f'deg-{dd} cone ' + ','.join(V[i] for i in tri), square=(dd <= 6))
if which in ('all', 'subsets'):
    run(core.parse('X0**6 + X1**6 + Y0**6 + Y1**6 + X0**3*X1*Y0*Y1/2'), 'edge', 'F|_{X2=Y2=0}')
    run(core.parse('X1**6 + X2**6 + Y1**6 + Y2**6 + X1**3*X2*Y1*Y2/2'), 'edge', 'F|_{X0=Y0=0}')
    run(core.parse('X0**6 + X2**6 + Y0**6 + Y2**6 + X0*X2**3*Y0*Y2/2'), 'edge', 'F|_{X1=Y1=0}')
    for c in ['X0**3*X1*Y0*Y1', 'X1**3*X2*Y1*Y2', 'X0*X2**3*Y0*Y2']:
        run(core.parse(c), 'single-coupling', f'{c} (support = F - coupling/2)')
    for i in [1, 2]:
        run(core.parse(f'X{i}**6 + Y{i}**6'), 'pair-fermat', f'X{i}^6 + Y{i}^6')
    run(core.parse('X0**6 + X1**6 + Y0**6 + Y1**6'), 'pair-fermat2', 'X0^6+X1^6+Y0^6+Y1^6 (edge minus coupling)')
    run(core.parse('X0**3*X1*Y0*Y1 + X1**3*X2*Y1*Y2'), 'two-couplings', 'c01 + c12')
    run(core.parse('X0**3*X1*Y0*Y1 + X0*X2**3*Y0*Y2'), 'two-couplings', 'c01 + c20')
    run(core.parse('X1**3*X2*Y1*Y2 + X0*X2**3*Y0*Y2'), 'two-couplings', 'c12 + c20')
    run(core.parse('X0**6 + Y0**6 + X1**6 + Y1**6 + X2**6 + Y2**6 + X0**3*X1*Y0*Y1/2'), 'F-minus-2coup', 'F - (c12+c20)/2')
    run(core.parse('X0**6 + Y0**6 + X1**6 + Y1**6 + X2**6 + Y2**6 + X1**3*X2*Y1*Y2/2'), 'F-minus-2coup', 'F - (c01+c20)/2')
    run(core.parse('X0**6 + Y0**6 + X1**6 + Y1**6 + X2**6 + Y2**6 + X0*X2**3*Y0*Y2/2'), 'F-minus-2coup', 'F - (c01+c12)/2')
if which in ('all', 'jacobian'):
    Fi = [core.pderiv(F, i) for i in range(6)]
    x = [core.parse(v) for v in V]
    for i in range(6):
        run(Fi[i], 'partial', f'F_{V[i]}', e1=False)
    for i in range(6):
        for j in range(6):
            run(core.mul(x[j], Fi[i]), 'x_j*F_i', f'{V[j]}*F_{V[i]}', e1=False)
    combos = {
     'sumY_jF_Yj': [0, 0, 0, 1, 1, 1], 'Grad1': [0, 1, 0, 0, 1, 0], 'Grad2': [0, 0, 1, 0, 0, 1],
     'X0F0-X1F1': [1, -1, 0, 0, 0, 0], 'X0F0-Y0F3': [1, 0, 0, -1, 0, 0], 'X0F0+X1F1-2X2F2': [1, 1, -2, 0, 0, 0],
     'X0F0+Y1F4': [1, 0, 0, 0, 1, 0], 'X0F0+Y2F5': [1, 0, 0, 0, 0, 1], 'X0F0+X1F1+Y0F3+Y1F4': [1, 1, 0, 1, 1, 0],
     'X0F0+2X1F1+3X2F2': [1, 2, 3, 0, 0, 0], 'alt': [1, -1, 1, -1, 1, -1]}
    random.seed(7)
    for k in range(4):
        combos[f'random-{k}'] = [random.randint(-5, 5) for _ in range(6)]
    for name, a in combos.items():
        D = {}
        for j in range(6):
            if a[j]: D = core.add(D, core.scale(core.mul(x[j], Fi[j]), Fraction(a[j])))
        if D: run(D, 'radial', name, e1=False)
    for i in range(6):
        for j in range(i, 6):
            run(core.mul(Fi[i], Fi[j]), 'F_i*F_j', f'F_{V[i]}*F_{V[j]}', square=False, e1=False)
print('DONE', flush=True)
