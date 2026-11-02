"""Like obsys.py but builds the J_{2d} reducers once and then processes several sub-families (Klein sectors, Y-parities,
slot subsets, or the full family) of the E1-admissible space at degree d.
Usage: python3 obsys2.py <d> <p> [--kernel f.json ...] [--threads T] [--fams full,sec00,sec01,sec10,sec11,par0,par1] [--no-run]
       [--msolve-timeout SEC]"""
import sys, time, json, os, subprocess, random, re
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
from common import *
args = sys.argv[1:]
d = int(args[0]); p = int(args[1])
kfiles = []; threads = 4; run = True; fams = ['full']; mto = 7200
i = 2
while i < len(args):
    if args[i] == '--kernel':
        i += 1
        while i < len(args) and not args[i].startswith('--'): kfiles.append(args[i]); i += 1
        continue
    if args[i] == '--threads': threads = int(args[i + 1]); i += 2; continue
    if args[i] == '--fams': fams = args[i + 1].split(','); i += 2; continue
    if args[i] == '--msolve-timeout': mto = int(args[i + 1]); i += 2; continue
    if args[i] == '--no-run': run = False; i += 1; continue
    i += 1
F, A, q, rec = load()
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
def chi_w(i, m):
    return ((m[5] - (1 if i == 5 else 0)) % 2, (m[3] + m[4] - (1 if i in (3, 4) else 0)) % 2)
# all unknowns of the full family
Uall = []
for i in range(NV):
    for m in monomials(d - 5):
        if i == 5 and m[5] >= 1: continue
        Uall.append(('w', i, m, chi_w(i, m)))
Ks = []; Ns = []
for kf in kfiles:
    Jf = json.load(open(kf))
    for s in Jf['basis']:
        K = reduce_mod(core.parse(s), p)
        N, info = e1_solve(F, A, K, p)
        assert N is not None, f'kernel class from {kf} not E1-admissible mod {p}: {info}'
        Ks.append(K); Ns.append(N); Uall.append(('K', len(Ks) - 1, None, klein(next(iter(K)))))
print(f'[obsys2 d={d} p={p}] full family: {len(Uall)} unknowns ({len(Ks)} kernel classes)', flush=True)
red = {par: Reducer(F, 2 * d, p, par, verbose=True) for par in (0, 1)}
def DN(u):
    if u[0] == 'w':
        _, i, m, _ = u
        D = shift(Fi[i], m, 1, p); N = [{} for _ in range(NV)]; N[i] = shift(Ap, m, 1, p)
        return D, N
    return Ks[u[1]], Ns[u[1]]
def nf_std(P):
    if not P: return {}
    R = red[core.parity_of(P)]
    M = R.nf_rows([P])
    return {R.mons[j]: int(M[0, j]) for j in R.standard if int(M[0, j])}
# self-tests
if d >= 6:
    rng = random.Random(1)
    l = {m: rng.randrange(1, p) for m in monomials(d - 6)}
    Dt = scale(mul(l, Fp, p), 6, p)
    Nt = [mul(Ap, mul(l, {tuple(1 if k == i else 0 for k in range(NV)): 1}, p), p) for i in range(NV)]
    assert nf_std(B(Dt, Nt, p)) == {}, 'self-test failed: ob(6lF) != 0'
    print('[obsys2] self-test ob(6lF)=0 mod J: OK', flush=True)
for k, (K, N) in enumerate(zip(Ks, Ns)):
    rng = random.Random(2 + k)
    N2 = [dict(N[i]) for i in range(NV)]
    for i in range(NV):
        for j in range(i + 1, NV):
            m = rng.choice(monomials(d - 4)); c = rng.randrange(1, p)
            N2[i] = add(N2[i], shift(Fi[j], m, c, p), p); N2[j] = add(N2[j], shift(Fi[i], m, p - c, p), p)
    nfa = nf_std(B(K, N, p)); nfb = nf_std(B(K, N2, p))
    assert nfa == nfb, 'self-test failed: ob depends on N'
    log({'exp': 'obsys-kernel-ob', 'd': d, 'p': p, 'k': k, 'ob_K_zero_modJ': nfa == {}, 'ob_K_nterms': len(nfa),
         'rigor': 'modular; ob(K) != 0 mod J mod p shows K itself (with any E1 solution reducing mod p) fails E2 mod J mod p'})
    print(f'[obsys2] kernel class {k}: ob independent of N OK; ob(K)==0 mod J: {nfa == {}} ({len(nfa)} terms)', flush=True)
def select(fam):
    if fam == 'full': return list(range(len(Uall)))
    if fam.startswith('sec'):
        sec = (int(fam[3]), int(fam[4])); return [a for a, u in enumerate(Uall) if u[3] == sec]
    if fam.startswith('par'):
        par = int(fam[3]); return [a for a, u in enumerate(Uall) if (u[3][0] + u[3][1]) % 2 == par]
    if fam.startswith('slots'):
        sl = [int(c) for c in fam[5:]]; return [a for a, u in enumerate(Uall) if u[0] == 'w' and u[1] in sl]
    raise ValueError(fam)
for fam in fams:
    t_start = time.time()
    idx = select(fam); U = [Uall[a] for a in idx]; nU = len(U)
    tag = f'obsys2_d{d}_p{p}_{fam}'
    Ds = []; Nsol = []
    for u in U:
        D, N = DN(u); Ds.append(D); Nsol.append(N)
    divN = []
    for N in Nsol:
        s = {}
        for i in range(NV): s = add(s, pderiv(N[i], i, p), p)
        divN.append(s)
    dD = [[pderiv(D, i, p) for i in range(NV)] for D in Ds]
    def Bab(a, b):
        out = mul(Ds[a], divN[b], p)
        for i in range(NV):
            if Nsol[b][i]: out = sub(out, mul(Nsol[b][i], dD[a][i], p), p)
        return out
    pairs = [(a, b) for a in range(nU) for b in range(a, nU)]
    print(f'[{tag}] {nU} unknowns, {len(pairs)} pairs', flush=True)
    eqs = {}
    CH = 3000
    for c0 in range(0, len(pairs), CH):
        chunk = pairs[c0:c0 + CH]
        polys = {0: [], 1: []}; who = {0: [], 1: []}
        for (a, b) in chunk:
            Q = Bab(a, b)
            if a != b: Q = add(Q, Bab(b, a), p)
            if not Q: continue
            par = core.parity_of(Q); polys[par].append(Q); who[par].append((a, b))
        for par in (0, 1):
            if not polys[par]: continue
            R = red[par]; M = R.nf_rows(polys[par])
            for r, (a, b) in enumerate(who[par]):
                for j in R.standard:
                    v = int(M[r, j])
                    if v: eqs.setdefault((par, j), []).append((v, a, b))
        print(f'[{tag}] pairs {c0 + len(chunk)}/{len(pairs)} ({time.time() - t_start:.0f}s) eqs {len(eqs)}', flush=True)
    neq = len(eqs)
    # rank of the quadric space (number of independent quadrics) mod p
    pidx = {}
    for key in eqs:
        for (v, a, b) in eqs[key]: pidx.setdefault((a, b), len(pidx))
    if neq:
        flat = [0] * (neq * len(pidx))
        for r, key in enumerate(sorted(eqs)):
            for (v, a, b) in eqs[key]: flat[r * len(pidx) + pidx[(a, b)]] = v
        qrank = nmod_mat(neq, len(pidx), flat, p).rank()
    else:
        qrank = 0
    msf = f'{PK}/{tag}.ms'
    with open(msf, 'w') as f:
        f.write(','.join(f't{a}' for a in range(nU)) + '\n' + f'{p}\n')
        lines = []
        for key in sorted(eqs):
            lines.append('+'.join(f'{v}*' + (f't{a}^2' if a == b else f't{a}*t{b}') for (v, a, b) in eqs[key]))
        f.write(',\n'.join(lines) + '\n')
    json.dump({'nU': nU, 'eqs': [[list(t) for t in eqs[key]] for key in sorted(eqs)],
               'unknowns': [[u[0], u[1], list(u[2]) if u[2] else None, list(u[3])] for u in U]}, open(f'{PK}/{tag}.eqs.json', 'w'))
    rec_ = {'exp': 'obsys2-build', 'd': d, 'p': p, 'family': fam, 'unknowns': nU, 'nK': sum(1 for u in U if u[0] == 'K'),
            'equations': neq, 'independent_quadrics': qrank, 'quadratic_monomials': nU * (nU + 1) // 2,
            'kernel_files': kfiles, 'seconds': round(time.time() - t_start, 1)}
    log(rec_)
    if not run or nU == 0:
        continue
    t0 = time.time(); outf = f'{PK}/{tag}.gb1'
    cmd = ['msolve', '-t', str(threads), '-g', '1', '-f', msf, '-o', outf]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=mto); rc = r.returncode; err = r.stderr[-500:]
    except subprocess.TimeoutExpired:
        rc = 'timeout'; err = ''
    ms_sec = time.time() - t0
    lead = open(outf).read() if os.path.exists(outf) and rc == 0 else ''
    body = lead[lead.find('['):lead.rfind(']')] if '[' in lead else ''
    mons_lead = [m.strip() for m in body.strip('[]').replace('\n', '').split(',') if m.strip()]
    pure = {}; bydeg = {}
    for m in mons_lead:
        dg = sum(int(e) for e in re.findall(r'\^(\d+)', m)); bydeg[dg] = bydeg.get(dg, 0) + 1
        if '*' not in m:
            mm = re.match(r't(\d+)\^(\d+)$', m)
            if mm: pure[int(mm.group(1))] = min(pure.get(int(mm.group(1)), 10 ** 9), int(mm.group(2)))
    m_primary = bool(mons_lead) and all(a in pure for a in range(nU))
    res = {'exp': 'obsys2-msolve', 'd': d, 'p': p, 'family': fam, 'unknowns': nU, 'equations': neq, 'independent_quadrics': qrank,
           'lead_monomials': len(mons_lead), 'lead_by_degree': bydeg, 'm_primary_mod_p': m_primary,
           'n_vars_without_pure_power': nU - len(pure), 'max_pure_power': max(pure.values()) if pure else None,
           'msolve_seconds': round(ms_sec, 1), 'msolve_rc': rc, 'threads': threads, 'stderr': err,
           'conclusion': (f'RIGOROUS: no certificate of degree {d} with D in the sub-family {fam} of the E1-admissible space (mod F-multiples); '
                          'm-primary mod p => m-primary over Q') if m_primary else
                         ('msolve did not finish' if rc != 0 else 'NOT m-primary mod p: nontrivial E2-mod-J solutions mod p; analyse')}
    log(res)
    json.dump({**res, 'pure_powers': pure, 'missing': [a for a in range(nU) if a not in pure][:100]}, open(f'{PK}/{tag}.result.json', 'w'), indent=1)
