import sys; sys.path.insert(0,'/tmp/hodge/results/one-reflection-skew')
import w2skew as W, numpy as np, pickle, time
kf = pickle.load(open('w2ker5_d7_r0s1.pkl','rb'))
kx = W.e1_kernel5_exact(7, 0, 1)
Kf = kf['Kbasis']
Kx = np.array([[float(v) for v in row] for row in kx['Kexact']])
Qx, _ = np.linalg.qr(Kx.T)
res = Kf - (Kf @ Qx) @ Qx.T
print('(0,1): float kernel rows outside exact span (norms):', np.round(np.linalg.norm(res, axis=1), 12), flush=True)
t0 = time.time()
try:
    kx00 = W.e1_kernel5_exact(7, 0, 0, max_primes=10, verbose=True)
    print('(0,0) reconstructed with', kx00['primes_used'], 'primes; max height', max(max(abs(v.numerator), v.denominator) for row in kx00['Kexact'] for v in row), round(time.time()-t0,1), 's', flush=True)
except Exception as ex:
    print('(0,0) fails with 10 primes:', ex, round(time.time()-t0,1), 's', flush=True)
kf00 = pickle.load(open('w2ker5_d7_r0s0.pkl','rb'))
print('(0,0) float: kerD', kf00['kerD'], 'sv_gap', kf00['sv_gap'], 'nonJ_float', kf00['nonJ_float'], 'nonJ exact', kf00['nonJ'], flush=True)
