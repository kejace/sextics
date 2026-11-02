"""Positive control for the wave-2 sector varpro: on accepted instances with a known explicit witness (d, D),
build the E1-admissible family of D's degree and character sector, check D lies in it and has residual ~0,
then run random restarts and record how often the optimizer finds a rel<1e-10 point and whether
reconstruction + exact confirmation recovers a witness.  Usage: python3 w2_control.py <job> <iid> [restarts] [procs]"""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['OMP_NUM_THREADS'] = '1'; os.environ['MKL_NUM_THREADS'] = '1'
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
import multiprocessing as mp

job = sys.argv[1]; iid = sys.argv[2]
NRES = int(sys.argv[3]) if len(sys.argv) > 3 else 64
PROCS = int(sys.argv[4]) if len(sys.argv) > 4 else 8
METHOD = sys.argv[5] if len(sys.argv) > 5 else 'lbfgs'
MAXIT = int(sys.argv[6]) if len(sys.argv) > 6 else 400
OUT = f'{LANE}/wave2-search/{job}.jsonl'
inst0 = Instance(iid)
d, Dw, Nw, Qw = core.load_witness(iid)
sub = symmetry.stabilizer(inst0.G, [inst0.A, Dw])
print(f'[control] {iid}: diagonal group order {len(inst0.G.elements)} (mod {inst0.n}); stabilizer of (A, D_witness) has {len(sub)} elements', flush=True)
inst = Instance(iid, gens=sub)
lab = label_of_poly(inst, Dw)
print(f'[control] witness d={d} label (under the {len(inst.Gsub.elements)}-element stabilizer) = {lab}', flush=True)
assert lab is not None
fam = e1_family(inst, d, lab, verbose=True)
dm = e1_family_dim_modp(inst, d, lab)
print('[control] family K =', fam['K'], 'mod-p dim', dm['dim'], 'J_d sector dim', jacobian_dim_sector(inst, d, lab), 'F-mult dim', fmult_dim_sector(inst, d, lab), flush=True)
svp = SectorVarPro(inst, d, lab, fam['basis'], verbose=True)
aw, prj = svp.coords_of(Dw)
relw = svp.rel(aw)
print(f'[control] witness projection residual {prj:.2e}; rel residual at witness {relw:.3e}; spurious? {svp.is_spurious(aw)}', flush=True)
log_line(OUT, {'kind': 'control', 'instance': iid, 'd': d, 'sector': lab, 'K': fam['K'], 'dimP': svp.dimP, 'rows': svp.rows, 'cols': svp.ncols,
               'nN': svp.nN, 'nQ': svp.nQ, 'full_nN': svp.full_nN, 'full_nQ': svp.full_nQ, 'witness_proj_res': prj, 'witness_rel': relw, 'method': METHOD, 'maxiter': MAXIT,
               'modp_dim': dm['dim'], 'engine': ENGINE})

def work(seed):
    t0 = time.time()
    if METHOD == 'lm':
        cands = svp.search_lm(n_restarts=8, seed=seed, max_nfev=MAXIT)
    elif METHOD == 'sparse':
        cands = svp.search(n_restarts=8, seed=seed, maxiter=MAXIT, init='sparse')
    elif METHOD == 'deflate':
        mins = svp.explore_minima(rounds=4, restarts=8, seed=seed, maxiter=MAXIT)
        cands = [{'rel': m['rel'], 'a': m['a'], 'spurious': False, 'iters': m['round']} for m in mins] or [{'rel': 1.0, 'a': svp.P[:, 0], 'spurious': False, 'iters': 0}]
    else:
        cands = svp.search(n_restarts=8, seed=seed, maxiter=MAXIT)
    res = []
    for c in cands:
        rec = {'rel': c['rel'], 'spurious': c['spurious'], 'iters': c['iters'], 'method': METHOD}
        if c['rel'] < 1e-8 and not c['spurious']:
            a1, r1 = svp.polish(c['a'])
            rec['polished_rel'] = r1
            Dfl = svp.D_float(a1)
            att = reconstruct_and_confirm(inst, Dfl, log_path=None, tag=f'control-{iid}-seed{seed}')
            rec['confirm'] = [(x['max_den'], x['status']) for x in att]
            rec['solved'] = any(x['status'] == 'solved' for x in att)
            if rec['solved']:
                rec['D'] = [x['D'] for x in att if x['status'] == 'solved'][0][:300]
        res.append(rec)
    return seed, res, time.time() - t0

if __name__ == '__main__':
    seeds = list(range(NRES // 8))
    with mp.get_context('fork').Pool(PROCS) as pool:
        for seed, res, dt in pool.imap_unordered(work, seeds):
            nfound = sum(1 for r in res if r['rel'] < 1e-10 and not r['spurious'])
            nsolved = sum(1 for r in res if r.get('solved'))
            print(f'[control] seed {seed}: {len(res)} restarts, {nfound} with rel<1e-10, {nsolved} exactly solved, best rel {min(r["rel"] for r in res):.2e} ({dt:.0f}s)', flush=True)
            log_line(OUT, {'kind': 'control-batch', 'instance': iid, 'd': d, 'sector': lab, 'seed': seed, 'restarts': len(res),
                           'n_rel_lt_1e-10': nfound, 'n_solved': nsolved, 'rels': [r['rel'] for r in res],
                           'solved_D': [r.get('D') for r in res if r.get('solved')][:1], 'seconds': round(dt, 1), 'method': METHOD, 'engine': ENGINE})
    print('[control] done', flush=True)
