"""Homogeneous quadratic obstruction system  ob(D) = [B_D(N_D)] = 0 in (S/J)_{2d}  on the whole E1-admissible family
   D = sum_k lambda_k K_k + sum_i w_i F_i   (w_i in S_{d-5}, modulo the multiples (F)_d of F, which do not change X cap V(D)),
where K_k are the non-Jacobian E1-kernel classes of degree d (mod p) and N_D = sum_k lambda_k N_k + A w  solves E1.
ob(D) is independent of the E1 solution N_D (syzygies of the regular sequence F_i are Koszul and land in J), and every
certificate has ob(D) = 0.  The system is written in msolve format (characteristic p) and solved for its leading ideal.
If the leading ideal mod p contains a pure power of every variable, the ideal is m-primary mod p, hence (Hilbert-function
semicontinuity: rank_p <= rank_Q of every Macaulay matrix of the same integer generators) m-primary over Q: no nonzero
E1-admissible D of degree d outside (F) satisfies E2 mod J, so no certificate of degree d exists.
Usage: python3 obsys.py <d> <p> [--kernel f.json ...] [--threads T] [--no-run]"""
import sys, time, json, os, subprocess, random
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
from common import *
args = sys.argv[1:]
d = int(args[0]); p = int(args[1])
kfiles = []; threads = 4; run = True
i = 2
while i < len(args):
    if args[i] == '--kernel':
        i += 1
        while i < len(args) and not args[i].startswith('--'):
            kfiles.append(args[i]); i += 1
        continue
    if args[i] == '--threads': threads = int(args[i + 1]); i += 2; continue
    if args[i] == '--no-run': run = False; i += 1; continue
    i += 1
F, A, q, rec = load()
Fp = reduce_mod(F, p); Ap = reduce_mod(A, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
t_start = time.time()
tag = f'obsys_d{d}_p{p}'
# ---- unknown basis
U = []   # entries: ('w', i, m) or ('K', k)
for i in range(NV):
    for m in monomials(d - 5):
        if i == 5 and m[5] >= 1:
            continue   # dropped: (5, l*Y2) spans a complement of {l*z : l in S_{d-6}}  (D = 6 l F)
        U.append(('w', i, m))
nW = len(U)
Ks = []; Ns = []
for kf in kfiles:
    J = json.load(open(kf))
    for s in J['basis']:
        K = reduce_mod(core.parse(s), p)
        N, info = e1_solve(F, A, K, p)
        assert N is not None, f'kernel class from {kf} is not E1-admissible mod {p}: {info}'
        Ks.append(K); Ns.append(N)
        U.append(('K', len(Ks) - 1))
nU = len(U)
print(f'[{tag}] unknowns: {nW} (J-part, complement of F-multiples) + {len(Ks)} kernel classes = {nU}', flush=True)
# ---- D_a, N_a
def DN(a):
    u = U[a]
    if u[0] == 'w':
        _, i, m = u
        D = shift(Fi[i], m, 1, p)
        N = [{} for _ in range(NV)]; N[i] = shift(Ap, m, 1, p)
        return D, N
    return Ks[u[1]], Ns[u[1]]
Ds = []; Nsol = []
for a in range(nU):
    D, N = DN(a); Ds.append(D); Nsol.append(N)
divN = []
for N in Nsol:
    s = {}
    for i in range(NV): s = add(s, pderiv(N[i], i, p), p)
    divN.append(s)
dD = [[pderiv(D, i, p) for i in range(NV)] for D in Ds]
def Bab(a, b):
    """B_{D_a}(N_b) = D_a div N_b - N_b(D_a)."""
    out = mul(Ds[a], divN[b], p)
    for i in range(NV):
        if Nsol[b][i]:
            out = sub(out, mul(Nsol[b][i], dD[a][i], p), p)
    return out
# ---- self-tests: ob(6 l F) = 0 and independence of N (Koszul) for a kernel class
red = {par: Reducer(F, 2 * d, p, par, verbose=True) for par in (0, 1)}
def nf_std(P):
    par = core.parity_of(P)
    if P == {}: return None, {}
    R = red[par]
    M = R.nf_rows([P])
    return par, {R.mons[j]: int(M[0, j]) for j in R.standard if int(M[0, j])}
if d >= 6:
    rng = random.Random(1)
    l = {m: rng.randrange(1, p) for m in monomials(d - 6)}
    Dt = scale(mul(l, Fp, p), 6, p)
    Nt = [mul(Ap, mul(l, {tuple(1 if k == i else 0 for k in range(NV)): 1}, p), p) for i in range(NV)]
    bt = B(Dt, Nt, p)
    par, nf = nf_std(bt)
    assert nf == {}, 'self-test failed: ob(6 l F) != 0 mod J'
    print(f'[{tag}] self-test ob(6lF)=0 mod J: OK', flush=True)
if Ks:
    K, N = Ks[0], Ns[0]
    rng = random.Random(2)
    Kos = [[{} for _ in range(NV)] for _ in range(NV)]
    for i in range(NV):
        for j in range(i + 1, NV):
            if d - 4 >= 0:
                m = rng.choice(monomials(d - 4)); Kos[i][j] = {m: rng.randrange(1, p)}; Kos[j][i] = scale(Kos[i][j], p - 1, p)
    N2 = [add(N[i], sum((mul(Kos[i][j], Fi[j], p) for j in range(NV)), {}) if False else N[i], p) for i in range(NV)]
    for i in range(NV):
        acc = dict(N[i])
        for j in range(NV):
            if Kos[i][j]: acc = add(acc, mul(Kos[i][j], Fi[j], p), p)
        N2[i] = acc
    chk = sub(sum_poly := {}, {}, p)
    e1a = {}; e1b = {}
    for i in range(NV):
        e1a = add(e1a, mul(N[i], Fi[i], p), p); e1b = add(e1b, mul(N2[i], Fi[i], p), p)
    assert e1a == e1b == mul(Ap, K, p), 'E1 check failed'
    _, nfa = nf_std(B(K, N, p)); _, nfb = nf_std(B(K, N2, p))
    assert nfa == nfb, 'self-test failed: ob depends on N'
    print(f'[{tag}] self-test ob independent of N (Koszul shift): OK; ob(K_0) == 0 mod J: {nfa == {}}', flush=True)
    log({'exp': 'obsys-kernel-ob', 'd': d, 'p': p, 'kernel_files': kfiles, 'ob_K0_zero_modJ': nfa == {}, 'ob_K0_nterms': len(nfa)})
# ---- all pairs
t0 = time.time()
pairs = [(a, b) for a in range(nU) for b in range(a, nU)]
print(f'[{tag}] building {len(pairs)} pair polynomials of degree {2*d}', flush=True)
eqs = {}   # (par, std-col-index) -> list of (coef, a, b)
CH = 4000
nnz_pairs = 0
for c0 in range(0, len(pairs), CH):
    chunk = pairs[c0:c0 + CH]
    polys = {0: [], 1: []}; who = {0: [], 1: []}
    for (a, b) in chunk:
        Q = Bab(a, b)
        if a != b: Q = add(Q, Bab(b, a), p)
        if not Q: continue
        par = core.parity_of(Q)
        assert par is not None
        polys[par].append(Q); who[par].append((a, b))
    for par in (0, 1):
        if not polys[par]: continue
        R = red[par]
        M = R.nf_rows(polys[par])
        for r, (a, b) in enumerate(who[par]):
            row_nz = False
            for j in R.standard:
                v = int(M[r, j])
                if v:
                    eqs.setdefault((par, j), []).append((v, a, b)); row_nz = True
            nnz_pairs += row_nz
    print(f'[{tag}] pairs {c0 + len(chunk)}/{len(pairs)} done ({time.time() - t0:.0f}s), eqs so far {len(eqs)}', flush=True)
neq = len(eqs)
print(f'[{tag}] {neq} nonzero quadrics in {nU} unknowns; pairs with nonzero ob: {nnz_pairs}', flush=True)
# ---- msolve file
msf = f'{PK}/{tag}.ms'
with open(msf, 'w') as f:
    f.write(','.join(f't{a}' for a in range(nU)) + '\n')
    f.write(f'{p}\n')
    lines = []
    for key in sorted(eqs):
        terms = []
        for (v, a, b) in eqs[key]:
            mon = f't{a}^2' if a == b else f't{a}*t{b}'
            terms.append(f'{v}*{mon}')
        lines.append('+'.join(terms))
    f.write(',\n'.join(lines) + '\n')
json.dump({'d': d, 'p': p, 'unknowns': [list(u) if u[0] == 'K' else [u[0], u[1], list(u[2])] for u in U], 'nW': nW, 'nK': len(Ks),
           'n_equations': neq, 'kernel_files': kfiles, 'build_seconds': round(time.time() - t_start, 1),
           'ring_dims_2d': {par: red[par].ring_dim for par in (0, 1)}}, open(f'{PK}/{tag}.meta.json', 'w'), indent=1)
log({'exp': 'obsys-build', 'd': d, 'p': p, 'unknowns': nU, 'nW': nW, 'nK': len(Ks), 'equations': neq, 'msolve_file': msf,
     'seconds': round(time.time() - t_start, 1)})
if not run:
    sys.exit(0)
# ---- msolve: leading ideal
t0 = time.time()
outf = f'{PK}/{tag}.gb1'
cmd = ['msolve', '-t', str(threads), '-g', '1', '-f', msf, '-o', outf]
print('[msolve]', ' '.join(cmd), flush=True)
r = subprocess.run(cmd, capture_output=True, text=True)
ms_sec = time.time() - t0
print(r.stdout[-2000:], r.stderr[-2000:], flush=True)
lead = open(outf).read() if os.path.exists(outf) else ''
# parse leading monomials
import re
body = lead[lead.find('['):lead.rfind(']')] if '[' in lead else ''
mons_lead = [m.strip() for m in body.strip('[]').replace('\n', '').split(',') if m.strip()]
pure = {}
for m in mons_lead:
    if '*' not in m:
        mm = re.match(r't(\d+)(?:\^(\d+))?$', m)
        if mm:
            a = int(mm.group(1)); e = int(mm.group(2) or 1)
            pure[a] = min(pure.get(a, 10**9), e)
m_primary = all(a in pure for a in range(nU))
missing = [a for a in range(nU) if a not in pure]
res = {'exp': 'obsys-msolve', 'd': d, 'p': p, 'unknowns': nU, 'equations': neq, 'lead_monomials': len(mons_lead),
       'm_primary_mod_p': m_primary, 'n_vars_without_pure_power': len(missing), 'max_pure_power': max(pure.values()) if pure else None,
       'msolve_seconds': round(ms_sec, 1), 'msolve_rc': r.returncode, 'threads': threads,
       'conclusion': ('RIGOROUS: no certificate with deg D = %d (E1 family mod (F) has only the zero solution of E2 mod J over C)' % d) if m_primary
                     else 'NOT m-primary mod p: nontrivial E2-mod-J solutions exist mod p; analyse'}
log(res)
json.dump({**res, 'missing_vars': missing[:50], 'pure_powers': pure}, open(f'{PK}/{tag}.result.json', 'w'), indent=1)
