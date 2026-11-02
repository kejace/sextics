"""JOB excl-dp-relax*: Hessian/E2-mod-J relaxation systems and msolve Groebner runs.

usage: python3 w2_relax.py <d> <sector> [--d8] [--threads T] [--timeout S] [--no-msolve]
  sector: e (Y-even D: eps=[0]), o (Y-odd D: eps=[1]), m (mixed: eps=[0,1]).
  --d8: include the D8 line (lam * D8 + w(F)); only meaningful for d = 8 (even / mixed).
Writes relax_d<d><sector>[_D8].ms (msolve input), .gb (msolve output), .json (receipt); logs exp='relax'.
Rigor: V(I) = {0} over the algebraic closure of GF(p)  ==>  V(I_Q) = {0}  (Hilbert-function argument; needs good
reduction of J in degree 2d, recorded).  Zero-dimensionality is read off the reduced Groebner basis: every
variable has a pure power among the leading monomials (grevlex)."""
import sys, os, time, json, subprocess, re
import numpy as np
import w2lib as W
from hodge import core

args = sys.argv[1:]
d = int(args[0]); sector = args[1]
with_d8 = '--d8' in args
threads = int(args[args.index('--threads') + 1]) if '--threads' in args else 16
timeout = int(args[args.index('--timeout') + 1]) if '--timeout' in args else 6 * 3600
run_msolve = '--no-msolve' not in args
p = int(args[args.index('--prime') + 1]) if '--prime' in args else W.P1
eps = {'e': [0], 'o': [1], 'm': [0, 1]}[sector]
tag = f'relax_d{d}{sector}' + ('_D8' if with_d8 else '') + (f'_p{p}' if p != W.P1 else '')
job = os.environ.get('W2JOB', tag)

D8 = None; u8 = None
if with_d8:
    D8 = core.parse(json.load(open(f'{W.LANE}/D8_e1kernel.json'))['D'])
    try:
        u8s = json.load(open(f'{W.PKT}/d8_u8_exact.json'))['u8']
        u8 = [core.reduce_mod(core.parse(s), p) for s in u8s]
    except Exception:
        u8 = None

t0 = time.time()
Rx = W.Relaxation(d, eps, p, with_D8=D8, u8=u8, euler_complement=True, verbose=True)
good = all(Rx.reds[par].ring_dim == W.koszul_hilbert(2 * d, par) for par in (0, 1))
E = Rx.independent_equations()
neq = E.shape[0]
npairs = len(Rx.pairs)
rec = {'d': d, 'sector': sector, 'with_D8': with_d8, 'prime': p, 'unknowns': Rx.n, 'pair_monomials': npairs,
       'raw_equations': Rx.neq_raw, 'independent_equations': int(neq), 'good_reduction_2d': good,
       'build_seconds': round(time.time() - t0, 1), 'names_head': Rx.names[:5]}
trivial = (neq == npairs)
rec['all_quadratic_monomials_in_I'] = trivial
print(rec, flush=True)
mspath = f'{W.PKT}/{tag}.ms'
Rx.write_msolve(mspath)
json.dump(rec, open(f'{W.PKT}/{tag}.json', 'w'), indent=1)
W.log('relax', family=tag, d=d, prime=p, result=('I contains m^2: V={0} (no msolve needed)' if trivial else 'system written'),
      rank=int(neq), seconds=rec['build_seconds'], note=f'unknowns {Rx.n}, pair monomials {npairs}, good reduction {good}')

if trivial or not run_msolve:
    print('DONE (trivial)' if trivial else 'DONE (no msolve)', flush=True)
    sys.exit(0)

# ------------------------------------------------------------------ msolve
gbpath = f'{W.PKT}/{tag}.gb'
cmd = ['msolve', '-f', mspath, '-o', gbpath, '-t', str(threads), '-g', '2', '-v', '1']
print('running', ' '.join(cmd), flush=True)
t1 = time.time()
try:
    pr = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    status = pr.returncode
    tail = (pr.stdout + pr.stderr)[-3000:]
except subprocess.TimeoutExpired:
    status = 'timeout'; tail = ''
secs = round(time.time() - t1, 1)
rec['msolve_status'] = status; rec['msolve_seconds'] = secs; rec['msolve_tail'] = tail
open(f'{W.PKT}/{tag}.msolve.log', 'w').write(tail)

verdict = 'msolve failed'
if status == 0 and os.path.exists(gbpath):
    txt = open(gbpath).read()
    # parse polynomials: msolve prints "[" poly "," poly ... "]"; strip brackets/whitespace
    body = '\n'.join(l for l in txt.splitlines() if not l.startswith('#')).strip()
    if body.startswith('['):
        body = body[1:]
    if body.endswith(']:') or body.endswith(']'):
        body = body.rstrip(':').rstrip(']')
    polys = [s.strip() for s in body.replace('\n', '').split(',') if s.strip()]
    names = Rx.names
    nidx = {nm: i for i, nm in enumerate(names)}
    n = len(names)

    def lead_monomial(poly):
        best = None
        for term in re.split(r'(?=[+-])', poly.replace(' ', '')):
            if not term or term in ('+', '-'):
                continue
            expo = [0] * n
            for fac in term.lstrip('+-').split('*'):
                if not fac or re.fullmatch(r'\d+', fac):
                    continue
                if '^' in fac:
                    v, e = fac.split('^'); expo[nidx[v]] += int(e)
                else:
                    expo[nidx[fac]] += 1
            key = (sum(expo), tuple(-x for x in reversed(expo)))   # grevlex: degree, then reverse-lex on reversed exponents
            if best is None or key > best[0]:
                best = (key, tuple(expo))
        return best[1] if best else None

    leads = [lead_monomial(q) for q in polys]
    leads = [l for l in leads if l is not None]
    pure = {}
    for l in leads:
        nz = [i for i, e in enumerate(l) if e]
        if len(nz) == 1:
            pure[nz[0]] = min(pure.get(nz[0], 10 ** 9), l[nz[0]])
    zero_dim = all(i in pure for i in range(n))
    rec['gb_size'] = len(polys); rec['gb_max_degree'] = max((sum(l) for l in leads), default=0)
    rec['pure_powers'] = len(pure); rec['zero_dimensional'] = zero_dim
    rec['gb_is_unit'] = (len(polys) == 1 and polys[0].strip() in ('1', '-1'))
    if zero_dim:
        verdict = 'V(I) = {0}: no E1-admissible D of this shape passes E2 mod J (rigorous over Q given good reduction)'
    else:
        # dimension of the leading-term ideal via a quick Singular call
        mono = ','.join('*'.join(f'{names[i]}^{e}' for i, e in enumerate(l) if e) for l in leads)
        sing = f'ring r = {p}, ({",".join(names)}), dp; ideal L = {mono}; int dm = dim(L); print(dm); quit;'
        try:
            sp = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=1800)
            dm = int(sp.stdout.strip().splitlines()[-1])
        except Exception as ex:
            dm = f'err {ex}'
        rec['dim_V'] = dm
        verdict = f'positive-dimensional: dim V(I) = {dm} (affine); needs analysis'
rec['verdict'] = verdict
json.dump(rec, open(f'{W.PKT}/{tag}.json', 'w'), indent=1)
W.log('relax_msolve', family=tag, d=d, prime=p, result=verdict, seconds=secs,
      note=f'status {status}; gb size {rec.get("gb_size")}; max lead degree {rec.get("gb_max_degree")}; zero_dim {rec.get("zero_dimensional")}')
print('DONE', verdict, flush=True)
