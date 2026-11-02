import pickle, time, json, sys, subprocess
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx
obj = pickle.load(open('d7_fam.pkl','rb'))
for k,(nm, Dd, N, den, chi, kind) in enumerate(obj['items']):
    if kind != 'extra' or k not in (90, 92, 122): continue
    out = w2x.test_fixed_D(ctx.from_dict(Dd), f'extra-dir d7 #{k} sector {tuple(chi)}', 'E1-kernel-nonJacobian-d7', want_exact_dual=True)
    print(k, chi, out['status'], out['rank'], out['rank_aug'], round(out['seconds'],1), flush=True)
print('extra tests done', flush=True)
subprocess.run([sys.executable, 'relax_bank.py', '7'])
