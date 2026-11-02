"""E2-mod-J quadratic relaxation over the E1-admissible family K_d + J(F)_d (gauge-fixed), solved by msolve/Singular.
Usage: python3 w2_relax.py <d> <parspec: 0|1|m> [--p P] [--threads T] [--solver msolve|singular|both] [--no-kernel]
       [--chart] [--tag name]
Reads kernel_d{d}_p{par}.json (rational basis preferred; else mod-p basis for the same prime).
Writes relax_d{d}_{parspec}.ms/.sing (homogeneous) and relax_d{d}_{parspec}_chart.ms (random affine chart), runs the
solver(s), stores raw outputs, logs a summary line."""
import sys, os, json, time, argparse, random, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *

def load_kernel(d, par, p):
    path = os.path.join(RES, f'kernel_d{d}_p{par}.json')
    if not os.path.exists(path):
        return None, 'missing'
    J = json.load(open(path))
    if J.get('basis_rational'):
        return [reduce_mod(core.parse(s), p) for s in J['basis_rational']], 'rational'
    if str(p) in J.get('basis_modp', {}):
        return [modp_poly_from_json(o) for o in J['basis_modp'][str(p)]], f'modp-{p}'
    return None, f'no basis for prime {p} (have {list(J.get("basis_modp", {}))})'

def gauge_indices(fam):
    idx = []
    pos = {(mem['j'], mem['m']): k for k, mem in enumerate(fam.members) if mem['tag'] == 'J'}
    for parD in fam.parities:
        for h in monomials(fam.d - 6, parD):
            key = (5, mon_add(h, (0, 0, 0, 0, 0, 1)))
            if key in pos:
                idx.append(pos[key])
    return sorted(set(idx))

def run_msolve(path, out, threads, timeout=None):
    t0 = time.time()
    try:
        r = subprocess.run(['msolve', '-t', str(threads), '-v', '1', '-f', path, '-o', out], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'rc': 'timeout', 'seconds': round(time.time() - t0, 1), 'out_head': 'TIMEOUT', 'stderr': '', 'stdout': '', 'out_len': 0}
    res = open(out).read() if os.path.exists(out) else ''
    return {'rc': r.returncode, 'stderr': r.stderr[-3000:], 'stdout': r.stdout[-3000:], 'out_head': res[:600],
            'out_len': len(res), 'seconds': round(time.time() - t0, 1)}

def run_singular(path, timeout):
    t0 = time.time()
    try:
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=timeout)
        return {'rc': r.returncode, 'stdout': r.stdout[-6000:], 'stderr': r.stderr[-2000:], 'seconds': round(time.time() - t0, 1)}
    except subprocess.TimeoutExpired:
        return {'rc': 'timeout', 'seconds': round(time.time() - t0, 1)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('parspec')
    ap.add_argument('--p', type=int, default=PGB[0]); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--solver', default='msolve'); ap.add_argument('--no-kernel', action='store_true')
    ap.add_argument('--tag', default=''); ap.add_argument('--singular-timeout', type=int, default=3600)
    ap.add_argument('--no-gauge', action='store_true'); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    set_threads(a.threads)
    F, A, q, rec = load()
    d, p = a.d, a.p
    pars = [0, 1] if a.parspec == 'm' else [int(a.parspec)]
    t0 = time.time()
    reps, nums, srcs = [], [], {}
    if not a.no_kernel:
        for par in pars:
            K, src = load_kernel(d, par, p)
            srcs[par] = src
            if K is None:
                print(f'[relax] kernel d={d} par={par}: {src}', file=sys.stderr)
                if src == 'missing' and d <= 6:
                    continue
                sys.exit(2)
            if K:
                N = particular_numerators(F, A, K, p)
                assert all(n is not None for n in N), 'kernel rep has no numerator (not E1-admissible mod p?)'
                reps += K; nums += N
    fam = Family(F, A, d, kernel_reps=reps, kernel_nums=nums, parities=pars)
    nK = len(reps)
    terms = relaxation_terms(fam, p)
    eqs, jrefs = reduce_terms(fam, terms, p)
    K = fam.K
    gidx = [] if a.no_gauge else gauge_indices(fam)
    keep = [k for k in range(K) if k not in set(gidx)]
    ren = {k: i for i, k in enumerate(keep)}
    eqs2 = []
    for eq in eqs:
        e2 = {}
        for (k, l), v in eq.items():
            if k in ren and l in ren:
                e2[(ren[k], ren[l])] = v
        if e2:
            eqs2.append(e2)
    n = len(keep)
    labels = [fam.members[k]['label'] for k in keep]
    base = os.path.join(RES, f'relax_d{d}_{a.parspec}' + (f'_{a.tag}' if a.tag else ''))
    write_msolve(base + '.ms', eqs2, n, p)
    write_singular(base + '.sing', eqs2, n, p, extra='if (dim(G) >= 1 and dim(G) <= 3) { LIB "primdec.lib"; list L = minAssGTZ(G); print("ncomp = " + string(size(L))); L; }\n')
    # random affine chart
    rng = random.Random(a.seed)
    chart = list(eqs2) + [{(): p - 1, **{(i,): rng.randrange(1, p) for i in range(n)}}]
    write_msolve(base + '_chart.ms', chart, n, p)
    write_singular(base + '_chart.sing', chart, n, p, extra='if (dim(G) == 0) { LIB "solve.lib"; print("vdim0 = " + string(vdim(G))); }\n')
    with open(base + '_meta.json', 'w') as f:
        json.dump({'d': d, 'parspec': a.parspec, 'p': p, 'K': K, 'nK': nK, 'gauge_removed': len(gidx), 'n': n,
                   'labels': labels, 'n_equations': len(eqs2), 'kernel_sources': srcs}, f)
    logrec = {'kind': 'relaxation', 'd': d, 'parspec': a.parspec, 'p': p, 'K_total': K, 'n_kernel': nK,
              'gauge_removed': len(gidx), 'n_vars': n, 'n_equations': len(eqs2), 'kernel_sources': srcs,
              't_build': round(time.time() - t0, 1), 'tag': a.tag}
    print(json.dumps(logrec), flush=True)
    if a.solver in ('msolve', 'both'):
        r1 = run_msolve(base + '_chart.ms', base + '_chart.msout', a.threads)
        logrec['msolve_chart'] = r1
        print('[msolve chart]', json.dumps(r1)[:800], flush=True)
        r0 = run_msolve(base + '.ms', base + '.msout', a.threads)
        logrec['msolve_hom'] = r0
        print('[msolve hom]', json.dumps(r0)[:800], flush=True)
    if a.solver in ('singular', 'both'):
        r2 = run_singular(base + '.sing', a.singular_timeout)
        logrec['singular_hom'] = r2
        print('[singular hom]', json.dumps(r2)[:1500], flush=True)
    logrec['seconds'] = round(time.time() - t0, 1)
    w2log(logrec)
    print(json.dumps(logrec, default=str)[:3000], flush=True)

if __name__ == '__main__':
    main()
