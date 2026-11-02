"""Control: on the pure Fermat sextic, which simple D serve X0^3 Y0^3? (menu item 4: deformation base point)"""
import sys, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
F0 = core.parse('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
A = core.parse('X0**3*Y0**3')
cands = ['X0**4 - X0**2*Y0**2 + Y0**4', 'X0**4 + Y0**4', 'X0**2*Y0**2', 'X0**2 + Y0**2', 'X0**2 - Y0**2', 'X0*Y0',
         'X0**4 + X1**4 + X2**4 + Y2**4', 'X0**4 + Y0**4 + Y2**4', 'X1**4 + X2**4 + Y1**4 + Y2**4',
         'X0**6 + Y0**6', 'X0**6 + Y0**6 + Y2**6', 'X1**6 + X2**6 + Y1**6 + Y2**6', 'X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6',
         'X0**6 + Y0**6 + Y1**6 + Y2**6', 'X1**6 + X2**6 + Y1**6', 'X1**6 + X2**6', 'X1**6 + Y1**6', 'X1**6 + 2*X2**6 + 3*Y1**6 + 5*Y2**6']
for e in cands:
    D = core.parse(e); t0 = time.time()
    ok, info = fixed_d.e1_solvable(F0, A, 1, D)
    st, data = fixed_d.solve_fixed_D(F0, A, 1, D, max_primes=3)
    print(f'{e:60s} d={core.degree(D)} e1={ok} status={st} rank={data.get("rank")} nullity={data.get("nullity")} ({time.time()-t0:.1f}s)', flush=True)
    if st == 'solved':
        print('   N terms', [len(n) for n in data['N']], 'Q terms', len(data['Q']))
