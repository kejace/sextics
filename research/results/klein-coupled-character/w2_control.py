"""Positive-control test of w2lib.KoszulVarPro on accepted d=6 witnesses (family = {D*} + J_6 generators).
Compares optimizers: L-BFGS (cutoff), varpro-LM (cutoff), varpro-LM (Tikhonov), joint LM.  Output: w2_control.out"""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character')
import w2lib as L
from hodge import core
import numpy as np
out = open('/tmp/hodge/results/klein-coupled-character/w2_control.out', 'a')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.write(s + '\n'); out.flush()
iid = sys.argv[1] if len(sys.argv) > 1 else 'hesse-triangle'
nres = int(sys.argv[2]) if len(sys.argv) > 2 else 10
variants = sys.argv[3].split(',') if len(sys.argv) > 3 else ['cutoff1e-9', 'tikh1e-6', 'tikh1e-8']
methods = sys.argv[4].split(',') if len(sys.argv) > 4 else ['lbfgs', 'lm', 'joint']
reduce_F = (sys.argv[5] != '0') if len(sys.argv) > 5 else True
F, A, q, rec = core.load_instance(iid)
d, Dw, Nw, Qw = core.load_witness(iid)
P('== control', iid, 'd', d, time.strftime('%H:%M:%S'))
t = time.time(); Nf, relf = L.float_e1_lift(F, A, Dw, restrict=False); P('float lift rel', relf, f'{time.time()-t:.1f}s')
fam = [(Dw, Nf, ('W',))] + L.jacobian_family(F, A, d, core.parity_of(Dw), None, restrict=False)
e0 = np.zeros(len(fam)); e0[0] = 1
VARS = {'cutoff1e-9': dict(cond=1e-9), 'cutoff1e-12': dict(cond=1e-12), 'tikh1e-6': dict(reg=1e-6), 'tikh1e-8': dict(reg=1e-8), 'tikh1e-10': dict(reg=1e-10)}
for name in variants:
    kw = VARS[name]
    kv = L.KoszulVarPro(F, A, fam, restrict=False, k=None, verbose=True, reduce_F=reduce_F, **kw)
    name = name + ('-red' if reduce_F else '-full')
    P(f'[{name}] rel at witness direction: {kv.rel(e0)[0]:.3e}')
    rng = np.random.default_rng(1)
    for sz in [0.3, 1.0]:
        a0 = e0 + sz * rng.standard_normal(kv.K)
        r = kv.search_lm(n_restarts=1, init=a0, max_nfev=200)[0]
        P(f'[{name}] varpro-LM from perturbation {sz}: rel={r["rel"]:.2e} nfev={r["iters"]} status={r["status"]}')
    for meth in methods:
        t = time.time()
        if meth == 'lbfgs': res = kv.search(n_restarts=nres, seed=int(os.environ.get('W2_SEED', '7')), maxiter=500)
        elif meth == 'lm': res = kv.search_lm(n_restarts=nres, seed=int(os.environ.get('W2_SEED', '7')), max_nfev=300)
        elif meth == 'clm': res = kv.search_clm(n_restarts=nres, seed=int(os.environ.get('W2_SEED', '7')), max_nfev=300)
        elif meth == 'two':
            res, stats = kv.search_two_stage(n_restarts=nres, seed=int(os.environ.get('W2_SEED', '7')), max_nfev=300, keep_frac=0.5)
            P(f'[{name}] two-stage relax stats: {stats}')
        else: res = kv.search_joint(n_restarts=nres, seed=int(os.environ.get('W2_SEED', '7')), max_nfev=300)
        hits = sum(1 for c in res if c['rel'] < 1e-14)
        rels = sorted(c['rel'] for c in res)
        P(f'[{name}] {meth}: hits {hits}/{len(res)} best {rels[0]:.2e} median {rels[len(rels)//2]:.2e} time {time.time()-t:.1f}s nfev/restart {np.mean([c["iters"] for c in res]):.0f}')
        if hits:
            best = res[0]
            cands = kv.reconstruct_D(best['a'], max_den=5000)
            Dr = L.rem_F_exact(Dw, F)
            def norm(Pp):
                m = max(Pp.values(), key=abs); return {e: v / m for e, v in Pp.items()}
            P(f'[{name}] {meth}: reconstruction candidates {len(cands)} match rem_F(D*): {any(norm(c)==norm(Dr) for c in cands)}; soldim {kv.solution_set_dim(best["a"])[0]}')
            if cands:
                st, data = L.confirm_exact(F, A, cands[0], restrict=False)
                P(f'[{name}] {meth}: confirm_exact -> {st} rank={data.get("rank")} verify={data.get("verify")}')
P('== done', time.strftime('%H:%M:%S'))
