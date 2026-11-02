import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
F, A, q, rec = core.load_instance('cyclic-triple-residual')
t0=time.time()
D = core.parse('1')
ok, info = fixed_d.e1_solvable(F, A, q, D)
print('D=1 e1_solvable:', ok, info, f'{time.time()-t0:.1f}s')
t0=time.time()
status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True, max_primes=1)
print('D=1 solve:', status, {k:v for k,v in data.items() if k not in ('N','Q')}, f'{time.time()-t0:.1f}s')
# timing on a degree-2 and degree-3 D
for Dexpr in ['X0**2 + X1**2', 'X0**3 + X1**3 + X2**3']:
    D = core.parse(Dexpr)
    t0=time.time()
    ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1=time.time()
    print(Dexpr, 'e1:', ok, info, f'{t1-t0:.1f}s')
    status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True, max_primes=1)
    print('   solve:', status, {k:v for k,v in data.items() if k not in ('N','Q')}, f'{time.time()-t1:.1f}s')
