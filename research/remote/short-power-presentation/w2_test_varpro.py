import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2_varpro import *
iid = sys.argv[1] if len(sys.argv) > 1 else 'fermat-pair-rational-d4-strict'
F, A, q, rec = core.load_instance(iid)
w = json.load(open(f'/home/kejace/devel/sextics/data/{iid}/witness.json'))
Dw = core.parse(w['D']); d = core.degree(Dw); parD = core.parity_of(Dw)
print('instance', iid, 'q', q, 'd', d, 'parD', parD, flush=True)
# exact verification of the witness with the engine
Nw = [core.parse(s) for s in w['N']]; Qw = core.parse(w['Q'])
print('verify_witness:', core.verify_witness(F, A, q, Dw, Nw, Qw) if hasattr(core, 'verify_witness') else 'n/a', flush=True)
t0 = time.time()
p1, p2 = PRIMES[0], PRIMES[1]
rb = []
for p in (p1, p2):
    Fp = core.reduce_mod(F, p); Ap = core.reduce_mod(A, p)
    reps, info, red = kernel_reps(Fp, Ap, d, parD, p, verbose=True)
    rb.append(reps)
rat, bad = crt_reconstruct(rb, [p1, p2])
print('kernel reps', len(rat), 'reconstructed ok' if rat else 'FAILED', f'{time.time()-t0:.0f}s', flush=True)
mems = J_members(F, A, d, parD) + kernel_members(F, A, d, parD, rat)
fam = ReducedVarPro(F, A, mems)
# coordinates of the witness D in the family
aw, res, rk, s = np.linalg.lstsq(fam.UD, fam.SD.vec(Dw), rcond=None)
print('witness in family: resid', np.linalg.norm(fam.UD @ aw - fam.SD.vec(Dw)), flush=True)
aw /= np.linalg.norm(aw)
ev = fam.evaluate(aw)
print('f at witness:', ev['f'], 'rel', ev['rel'], 'b2', ev['b2'], flush=True)
# consistency: CK c == T_D(c)
D = fam.UD @ aw; Di = [fam.UDi[i] @ aw for i in range(6)]
CK = fam.build_CK(D, Di)
c = np.random.default_rng(0).standard_normal(fam.nK)
V = fam.V_from_c(c); divV = sum(fam.Dmaps_N[i] @ V[i] for i in range(6))
T = fam.pm1.prod(D, divV) - sum(fam.pm2[i].prod(V[i], Di[i]) for i in range(6))
print('CK c vs T_D(c):', np.linalg.norm(CK @ c - T) / np.linalg.norm(T), flush=True)
# gradient check
u0 = np.random.default_rng(1).standard_normal(fam.dimP)
f0, g0 = fam.value_and_grad(u0)
h = 1e-6; du = np.random.default_rng(2).standard_normal(fam.dimP)
f1, _ = fam.value_and_grad(u0 + h * du); f2, _ = fam.value_and_grad(u0 - h * du)
print('grad check:', (f1 - f2) / (2 * h), 'vs', g0 @ du, flush=True)
# random point after projection: is witness in P-range?
print('witness component outside P:', np.linalg.norm(aw - fam.P @ (fam.P.T @ aw)), flush=True)
cands = fam.search(n_restarts=int(sys.argv[2]) if len(sys.argv) > 2 else 20, seed=0, verbose=True)
print('best rel', cands[0]['rel'], 'spur', cands[0]['spurious'])
best = np.array(cands[0]['a'])
for fr in fam.reconstruct(best, max_den=1000):
    Dx = fam.D_exact(fr); print('reconstructed D:', core.to_str(Dx)[:200])
    st, data = fam.confirm(fr); print('confirm:', st)
    break
