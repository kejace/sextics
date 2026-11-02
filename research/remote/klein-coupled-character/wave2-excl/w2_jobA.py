"""Job A (remote): (1) engine validation on accepted witnesses, (2) colon ideal / E1-kernel sweep d=1..12 both parities,
all characters, prime 1000003 (re-verification of dossier §1.7 / attack note §2.1 with a new prime and new code),
(3) Hessian-relaxation systems: d=5,6 exact Macaulay proof; d=7,8 eigen sectors solved with msolve in-job;
full (non-eigen) d=7,8 and mixed-parity systems written to msolve files for job B.
Outputs: jobA_*.json receipts, jobA_msolve/*.ms inputs, log lines in wave2-excl.log.jsonl."""
import sys, os, json, time, subprocess
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
P = int(sys.argv[1]) if len(sys.argv) > 1 else 1000003
F, A, q = load()
ctx = Ctx(F, A, P)
OUT = PACKET + f'/jobA_p{P}'
os.makedirs(OUT, exist_ok=True)
def jdump(obj, name):
    json.dump(obj, open(f'{OUT}/{name}', 'w'), indent=1, default=str)
def say(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)

# (1) engine validation: restricted solver must find accepted witnesses consistent and a perturbed D inconsistent
val = []
for iid in ['collision-h3', 'hesse-triangle', 'contact46006']:
    try:
        Fv, Av, qv, rec = core.load_instance(iid)
        d, Dv, Nv, Qv = core.load_witness(iid)
        t = time.time()
        st, data = symmetry.solve_fixed_D_symmetric(Fv, Av, qv, Dv, primes=[P, 1073741827], max_primes=2)
        D2 = dict(Dv); e0 = next(iter(D2)); D2[e0] = D2[e0] + 1
        st2, data2 = symmetry.solve_fixed_D_symmetric(Fv, Av, qv, D2, primes=[P], max_primes=1)
        val.append({'instance': iid, 'd': d, 'status_witness_D': st, 'rank': data.get('rank'), 'subgroup': data.get('subgroup_size'),
                    'status_perturbed_D': st2, 'seconds': round(time.time() - t, 1)})
        say('validation', val[-1])
    except Exception as ex:
        val.append({'instance': iid, 'error': repr(ex)}); say('validation error', iid, ex)
jdump(val, 'validation.json')
log({'exp': 'engine-validation', 'prime': P, 'result': val, 'job': 'jobA'})

# (2) colon / kernel sweep
sweep = {}
kernel_store = {}
for d in range(1, 13):
    t = time.time()
    row = {}
    for par in (0, 1):
        for k in range(7):
            cs = colon_sector(ctx, d, par, k, want_reps=(d >= 9))
            hf = hilbert_R_sector(d, par, k)
            row[f'{par},{k}'] = {'dimS': cs['dimS'], 'dimJ': cs['dimJ'], 'HF_R': hf, 'HF_ok': cs['dimS'] - cs['dimJ'] == hf,
                                 'dimColon': cs['dimColon'], 'kernel': cs['kernel']}
            if d >= 9 and cs['kernel']:
                kernel_store[f'{d},{par},{k}'] = [{str(m): v for m, v in P_.items()} for P_ in cs['reps']]
    tot = {par: sum(row[f'{par},{k}']['kernel'] for k in range(7)) for par in (0, 1)}
    hf_ok = all(v['HF_ok'] for v in row.values())
    sweep[d] = row
    say(f'd={d}: kernel even {[row[f"0,{k}"]["kernel"] for k in range(7)]} (tot {tot[0]}), odd {[row[f"1,{k}"]["kernel"] for k in range(7)]} (tot {tot[1]}), HF_ok={hf_ok} [{time.time()-t:.0f}s]')
    log({'exp': 'colon-kernel-sweep', 'd': d, 'prime': P, 'result': {'kernel_even': [row[f'0,{k}']['kernel'] for k in range(7)],
         'kernel_odd': [row[f'1,{k}']['kernel'] for k in range(7)], 'HF_ok': hf_ok,
         'dimJ_even': [row[f'0,{k}']['dimJ'] for k in range(7)], 'dimJ_odd': [row[f'1,{k}']['dimJ'] for k in range(7)]},
         'seconds': round(time.time() - t, 1), 'job': 'jobA',
         'note': 'dim (J:A)_d/J_d per (parity, sigma-character) mod p; HF_ok: dim J_d mod p equals the complete-intersection value in every sector (flatness at p)'})
    jdump(sweep, 'colon_sweep.json'); jdump(kernel_store, 'kernel_reps.json')

# engine cross-check of parity totals at d=9,10 (independent code path: hodge.jacobian.mult_kernel)
for d in (9, 10):
    for par in (0, 1):
        t = time.time()
        mk = jacobian.mult_kernel(F, A, d, P, parity=par)
        say(f'engine mult_kernel d={d} par={par}: ker[A]={mk["dim"]} J_d={mk["jacobian_dim"]} R_d={mk["ring_dim"]} [{time.time()-t:.0f}s]')
        log({'exp': 'engine-mult-kernel', 'd': d, 'par': par, 'prime': P, 'result': {'ker': mk['dim'], 'dimJ': mk['jacobian_dim'], 'dimR': mk['ring_dim']},
             'seconds': round(time.time() - t, 1), 'job': 'jobA', 'note': 'hodge.jacobian.mult_kernel (parity only), cross-check of the sector sums'})

# (3) Hessian relaxation
MS = PACKET + '/jobA_msolve'; os.makedirs(MS, exist_ok=True)
relax = {}
def run_msolve(path, threads=4, timeout=3600):
    out = path[:-3] + '.lead'
    t = time.time()
    try:
        r = subprocess.run(['msolve', '-f', path, '-o', out, '-g', '1', '-t', str(threads)], capture_output=True, text=True, timeout=timeout)
        rc = r.returncode
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'seconds': round(time.time() - t, 1)}
    if rc != 0 or not os.path.exists(out):
        return {'status': 'error', 'rc': rc, 'stderr': r.stderr[-300:], 'seconds': round(time.time() - t, 1)}
    return {'status': 'done', 'lead': out, 'seconds': round(time.time() - t, 1)}

def record(tag, s, extra):
    gauge_tested, gauge_bad = hess_check_gauge(ctx, s['d'], s['Dpars'], s['Dchars']) if s['d'] >= 6 else (0, 0)
    rec = {'tag': tag, 'd': s['d'], 'Dpars': s['Dpars'], 'Dchars': s['Dchars'], 'n_vars': s['n'], 'n_before_gauge': s['n_unknowns_before_gauge'],
           'n_gauge': s['n_gauge'], 'n_eqs': s['n_eqs'], 'targets': s['target_info'], 'gauge_dirs_tested': gauge_tested, 'gauge_violations': gauge_bad, 'prime': P}
    rec.update(extra)
    relax[tag] = rec
    jdump(relax, 'relaxation.json')
    log({'exp': 'hessian-relaxation', 'family': tag, 'd': s['d'], 'prime': P, 'result': rec.get('verdict'), 'rank': rec.get('macaulay', {}).get('rank'),
         'seconds': rec.get('seconds'), 'job': 'jobA', 'note': json.dumps({k: v for k, v in rec.items() if k not in ('targets',)}, default=str)[:900]})

# d = 5, 6: exact Macaulay-matrix proofs (degree 2 resp. 3)
for d in (5, 6):
    for pars, ptag in (([0], 'even'), ([1], 'odd'), ([0, 1], 'mixed')):
        t = time.time()
        s = hess_system(ctx, d, pars, range(7))
        mr = macaulay_rank(s, 2 if d == 5 else 3)
        verdict = 'EMPTY (all forms of degree %d in the ideal; projective variety empty over Qbar)' % mr['deg'] if mr['full_column_rank'] else 'not decided by Macaulay rank'
        say(f'relax d={d} {ptag}: n={s["n"]} eqs={s["n_eqs"]} macaulay={mr} -> {verdict} [{time.time()-t:.0f}s]')
        record(f'd{d}-{ptag}-all-chars', s, {'macaulay': mr, 'verdict': verdict, 'seconds': round(time.time() - t, 1)})
        write_msolve(s, f'{MS}/d{d}_{ptag}_all.ms')

# full (non-eigen) systems for job B
for d in (7, 8):
    for pars, ptag in (([0], 'even'), ([1], 'odd'), ([0, 1], 'mixed')):
        t = time.time()
        s = hess_system(ctx, d, pars, range(7))
        path = f'{MS}/d{d}_{ptag}_all.ms'
        write_msolve(s, path)
        say(f'wrote full system d={d} {ptag}: n={s["n"]} eqs={s["n_eqs"]} -> {path} [{time.time()-t:.0f}s]')
        record(f'd{d}-{ptag}-all-chars', s, {'msolve': 'deferred to job B', 'verdict': 'pending', 'seconds': round(time.time() - t, 1), 'file': path})
# d = 7, 8: eigen sectors (each character separately; pure parities and mixed parity), msolve in-job
for d in (7, 8):
    for k in range(7):
        for pars, ptag in (([0], 'even'), ([1], 'odd'), ([0, 1], 'mixed')):
            t = time.time()
            s = hess_system(ctx, d, pars, [k])
            path = f'{MS}/d{d}_{ptag}_k{k}.ms'
            write_msolve(s, path)
            res = run_msolve(path, threads=8, timeout=7200)
            verdict = None; lead = None
            if res['status'] == 'done':
                mons = parse_msolve_leading(res['lead'])
                irr, info = leading_ideal_irrelevant(mons, s['n'])
                cnt = standard_monomial_count(mons, s['n']) if irr and s['n'] <= 30 else None
                verdict = 'EMPTY (leading ideal irrelevant: zero-dimensional cone => projective variety empty over Qbar)' if irr else 'NONEMPTY mod p (positive-dimensional cone)'
                lead = {'n_lead': len(mons), 'irrelevant': irr, 'info': info if not irr else {'pure_powers': info['pure_powers']}, 'std_monomials': cnt}
            say(f'relax d={d} k={k} {ptag}: n={s["n"]} eqs={s["n_eqs"]} msolve={res["status"]} {res.get("seconds")}s -> {verdict}')
            record(f'd{d}-{ptag}-k{k}', s, {'msolve': res, 'lead': lead, 'verdict': verdict, 'seconds': round(time.time() - t, 1)})

say('JOBA_DONE')
