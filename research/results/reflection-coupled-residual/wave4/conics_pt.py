"""M5 step 1: weighted conics (e=2) on W = X/r in P(1^5,2) through a random point of W, over F_p.
W: Ft(x,u) = F(x, Y2) with u = Y2^2, x = (X0,X1,X2,Y0,Y1).  A weighted conic is tau -> (a(tau), b(tau)),
a_i quadratic, b quartic, with Ft(a,b) == 0 (degree-12 binary form, 13 coefficients).
Through a point (x0,u0) in W put the point at tau=(1:0): a_i = x0_i t0^2 + m_i t0 t1 + n_i t1^2,
b = u0 t0^4 + b1 t0^3 t1 + b2 t0^2 t1^2 + b3 t0 t1^3 + b4 t1^4; residual group (t0,t1)->(t0+e t1, c t1)
normalised by b1 = 0, m_0 = 1 (after a random linear change of the 5 weight-1 coordinates).
Unknowns m1..m4, n0..n4, b2,b3,b4 (12); equations: coefficients of t0^(12-k) t1^k, k=1..12 (12).
Usage: python3 conics_pt.py <name> <prime> <seed> [threads] [mode]   mode: pt (default) | hyp
  hyp: conics inside the random hyperplane section {x_4' = 0}, PGL2-normalised a_0 = t0 t1, a_1 = t0^2 + x t0 t1 + y t1^2
       (13 unknowns: x,y, a_2,a_3 (3 each), b (5); 13 equations).  Each conic counted twice (swap of the two roots).
Writes <name>.ms, <name>.param, <name>.json
"""
import sys, json, time, random, subprocess, re, os, shutil
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
import flint

name, p, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
threads = sys.argv[4] if len(sys.argv) > 4 else '4'
mode = sys.argv[5] if len(sys.argv) > 5 else 'pt'
random.seed(seed)
t0 = time.time()
F, A, q, rec = core.load_instance('reflection-coupled-residual')
# F -> Ft(x0..x4, u): Y2 exponent must be even
Ftd = {}
for m, cf in F.items():
    assert m[5] % 2 == 0
    key = (m[0], m[1], m[2], m[3], m[4], m[5] // 2)
    Ftd[key] = (cf.numerator * pow(cf.denominator, -1, p)) % p
XV = ['x0', 'x1', 'x2', 'x3', 'x4', 'u']
ctxW = flint.nmod_mpoly_ctx.get(tuple(XV), ordering=flint.Ordering.lex, modulus=p)
Ft = ctxW.from_dict(Ftd)
gW = ctxW.gens()
# random linear change on the 5 weight-1 coordinates
while True:
    G = [[random.randrange(p) for _ in range(5)] for _ in range(5)]
    if flint.nmod_mat(G, p).rank() == 5:
        break
FtG = Ft.compose(*[sum((G[a][b] * gW[b] for b in range(5)), ctxW.from_dict({})) for a in range(5)] + [gW[5]])

def cubic_root(x0):
    # u^3 + L(x0) u + M(x0) = 0 : evaluate FtG at (x0, u) as poly in u, find roots in F_p
    pu = FtG.compose(*[ctxW.from_dict({(0,)*6: int(v)}) for v in x0] + [gW[5]])
    d = pu.to_dict()
    coeffs = [0, 0, 0, 0]
    for m, cf in d.items():
        coeffs[m[5]] = int(cf)
    Pu = flint.nmod_poly(coeffs, p)
    rts = [int(r) for r, mult in Pu.roots()]
    return [r for r in rts if r != 0]

if mode == 'pt':
    while True:
        x0 = [random.randrange(1, p) for _ in range(5)]
        rts = cubic_root(x0)
        if rts:
            u0 = random.choice(rts)
            break
    UN = [f'm{i}' for i in range(1, 5)] + [f'n{i}' for i in range(5)] + ['b2', 'b3', 'b4']
    ALL = ['t0', 't1'] + UN
    ctxA = flint.nmod_mpoly_ctx.get(tuple(ALL), ordering=flint.Ordering.lex, modulus=p)
    gA = ctxA.gens()
    T0, T1 = gA[0], gA[1]
    def un(nm): return gA[ALL.index(nm)]
    const = lambda c: ctxA.from_dict({(0,) * len(ALL): int(c) % p})
    a = []
    for i in range(5):
        mi = const(1) if i == 0 else un(f'm{i}')
        a.append(const(x0[i]) * T0 ** 2 + mi * T0 * T1 + un(f'n{i}') * T1 ** 2)
    b = const(u0) * T0 ** 4 + un('b2') * T0 ** 2 * T1 ** 2 + un('b3') * T0 * T1 ** 3 + un('b4') * T1 ** 4
    info = {'x0': x0, 'u0': u0}
elif mode == 'hyp':
    UN = ['x', 'y'] + [f'a{i}_{j}' for i in (2, 3) for j in range(3)] + [f'b{j}' for j in range(5)]
    ALL = ['t0', 't1'] + UN
    ctxA = flint.nmod_mpoly_ctx.get(tuple(ALL), ordering=flint.Ordering.lex, modulus=p)
    gA = ctxA.gens()
    T0, T1 = gA[0], gA[1]
    def un(nm): return gA[ALL.index(nm)]
    const = lambda c: ctxA.from_dict({(0,) * len(ALL): int(c) % p})
    a = [T0 * T1, T0 ** 2 + un('x') * T0 * T1 + un('y') * T1 ** 2]
    for i in (2, 3):
        a.append(un(f'a{i}_0') * T0 ** 2 + un(f'a{i}_1') * T0 * T1 + un(f'a{i}_2') * T1 ** 2)
    a.append(ctxA.from_dict({}))  # x4' = 0 : the hyperplane section
    b = sum((un(f'b{j}') * T0 ** (4 - j) * T1 ** j for j in range(5)), ctxA.from_dict({}))
    info = {}
else:
    raise SystemExit('mode?')
# lift FtG to ctxA and compose
FtA = ctxA.from_dict({m + (0,) * (len(ALL) - 6): cf for m, cf in
                      {(mm[0], mm[1], mm[2], mm[3], mm[4], mm[5]): c for mm, c in FtG.to_dict().items()}.items()}) if False else None
# simpler: compose FtG directly with polynomials in ctxA by rebuilding term by term
comp = ctxA.from_dict({})
for m, cf in FtG.to_dict().items():
    term = const(cf)
    for i in range(5):
        if m[i]:
            term = term * a[i] ** m[i]
    if m[5]:
        term = term * b ** m[5]
    comp = comp + term
eqs = {}
for m, cf in comp.to_dict().items():
    k = m[1]  # exponent of t1
    rest = m[2:]
    eqs.setdefault(k, {})[rest] = cf
zero = ctxA.from_dict({})
eqlist = []
for k in sorted(eqs):
    e = ctxA.from_dict({(0, 0) + r: cf for r, cf in eqs[k].items()})
    if e != zero:
        eqlist.append(e)
print(f'[{name}] mode={mode} p={p} unknowns={len(UN)} equations={len(eqlist)} (t1-degrees {sorted(eqs)}) built in {time.time()-t0:.1f}s', flush=True)
def ms_str(e):
    s = []
    for m, cf in e.to_dict().items():
        term = [str(cf)]
        for k, ex in enumerate(m):
            if ex:
                term.append(f'{ALL[k]}^{ex}' if ex > 1 else ALL[k])
        s.append('*'.join(term))
    return '+'.join(s) if s else '0'
with open(f'{name}.ms', 'w') as f:
    f.write(','.join(UN) + '\n' + str(p) + '\n')
    f.write(',\n'.join(ms_str(e) for e in eqlist) + '\n')
summary = {'name': name, 'mode': mode, 'p': p, 'seed': seed, 'G': G, 'unknowns': UN, 'n_eqs': len(eqlist), **info}
json.dump(summary, open(f'{name}.json', 'w'))
if shutil.which('msolve') is None:
    print('msolve not found'); sys.exit(0)
t1 = time.time()
cmd = ['msolve', '-v', '1', '-t', threads, '-P', '1', '-f', f'{name}.ms', '-o', f'{name}.param']
r = subprocess.run(cmd, capture_output=True, text=True)
out = open(f'{name}.param').read() if r.returncode == 0 and os.path.exists(f'{name}.param') else ''
summary['msolve_rc'] = r.returncode; summary['msolve_seconds'] = time.time() - t1
summary['stderr_tail'] = r.stderr[-600:]
m = re.match(r'\s*\[\s*(-?\d+)\s*,', out)
if m:
    summary['dim'] = int(m.group(1))
    m2 = re.match(r'\s*\[\s*0\s*,\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,', out)
    if m2:
        summary['degree'] = int(m2.group(3))
summary['param_head'] = out[:160].replace('\n', ' ')
print(json.dumps({k: v for k, v in summary.items() if k != 'G'}), flush=True)
json.dump(summary, open(f'{name}.json', 'w'))
