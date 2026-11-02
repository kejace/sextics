"""Job: full E1 kernel dims at d=7 (both parities) for (a) other targets A' on F, (b) other coupling constants c, (c) the family F_s = eX+eY+B(sL,M).
Usage: python3 w2c_k7var.py <out.jsonl>"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2c_lib import *
import flint
try: flint.ctx.threads = int(os.environ.get('FLINT_THREADS', '8'))
except Exception: pass
out_path = sys.argv[1]
F, A, q, rec = core.load_instance(IID)
p = core.PRIMES[0]
def Fof(cc, s=Fraction(1)):
    sL = core.scale(L, s)
    U = core.add(sL, M); V = core.sub(sL, M)
    return core.add(core.add(eX, eY), core.scale(core.add(ppow(U, 6), ppow(V, 6)), cc))
assert core.sub(Fof(c), F) == {}
tasks = []
for name, Ap in [('X0^3Y1^3', 'X0**3*Y1**3'), ('X1^3Y1^3', 'X1**3*Y1**3'), ('X0X1X2Y0Y1Y2', 'X0*X1*X2*Y0*Y1*Y2'),
                 ('X0^4X1Y0', 'X0**4*X1*Y0'), ('L^3M^3', '(X0+2*X1+3*X2)**3*(Y0+3*Y1+5*Y2)**3'),
                 ('X0^3M^3', 'X0**3*(Y0+3*Y1+5*Y2)**3'), ('L^3Y0^3', '(X0+2*X1+3*X2)**3*Y0**3')]:
    tasks.append(('target', name, F, core.parse(Ap)))
for cc in [Fraction(1), Fraction(-1, 2), Fraction(1, 3), Fraction(5, 7)]:
    tasks.append(('coupling-c', str(cc), Fof(cc), A))
for s in [Fraction(2), Fraction(1, 2), Fraction(3, 5)]:
    tasks.append(('family-s', str(s), Fof(c, s), A))
for d in (7,):
    for kind, name, FF, AA in tasks:
        for par in (0, 1):
            t0 = time.time()
            out = jacobian.mult_kernel(FF, AA, d, p, parity=par, verbose=True)
            rec_out = {'exp': 'E1-kernel-variation', 'kind': kind, 'name': name, 'd': d, 'parD': par, 'prime': p,
                       'kernel_mod_J': out['dim'], 'kernel_dim': out['kernel_dim'], 'J_d': out['jacobian_dim'], 'R_d': out['ring_dim'],
                       'seconds': round(time.time() - t0, 1), 'job': 'construct-spp-k7var', 'result': 'kerJ=%d' % out['dim'],
                       'basis_nterms': [len(b) for b in out['basis']]}
            # structure test of basis elements: depends on Y only via M? on X only via L?
            st = []
            for b in out['basis']:
                dY = [core.pderiv(b, i, p) for i in range(3, 6)]
                viaM = core.sub(dY[1], core.scale(dY[0], 3, p), p) == {} and core.sub(dY[2], core.scale(dY[0], 5, p), p) == {}
                dX = [core.pderiv(b, i, p) for i in range(3)]
                viaL = core.sub(dX[1], core.scale(dX[0], 2, p), p) == {} and core.sub(dX[2], core.scale(dX[0], 3, p), p) == {}
                st.append({'Y_via_M': viaM, 'X_via_L': viaL})
            rec_out['structure'] = st
            with open(out_path, 'a') as f: f.write(json.dumps(rec_out) + '\n')
            print('DONE', kind, name, d, par, 'kerJ', out['dim'], st, flush=True)
print('ALL DONE', flush=True)
