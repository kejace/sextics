"""Wave-2 (excl lens): quadratic relaxation of E2 for ALL denominators of degree d with a given parity structure.

Usage: python3 excl_relax.py <d> <pars> [--pidx k] [--threads T] [--no-singular] [--no-msolve] [--no-slice] [--tag TAG]
   pars: e (Y-even D), o (Y-odd D), m (mixed D = D^+ + D^-)

Pipeline (all mod p = PRIMES30[pidx]):
  1. K_d^{par}/J_d^{par} basis (excl_lib.e1_kernel) and particular numerators N_a (N_a(F) = A D_a);
  2. the quadratic system Xi(lam, w) in J_{2d} (excl_lib.build_relaxation) in the Euler-slice coordinates;
  3. msolve -g 1 (leading ideal, F4, threads) -> Hilbert function of S/I_p by counting standard monomials;
     Singular std/dim/vdim/hilbPoly as an independent check;
  4. verdict 'proj-empty' iff the Hilbert function vanishes in some degree: then, by the Hilbert-function
     inequality h_{S/I_p} >= h_{S/I_Q}, no denominator of degree d with this parity structure exists over
     Q-bar except multiples of F (not coprime).
"""
import sys, time, json, os, subprocess, re, itertools, math
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
import excl_lib as W
from hodge import core


def parse_msolve_leading(path):
    txt = ''.join(l for l in open(path).read().splitlines() if not l.startswith('#'))
    txt = txt.replace(' ', '')
    txt = txt.strip().rstrip(':').strip()
    if txt.startswith('['):
        txt = txt[1:]
    if txt.endswith(']'):
        txt = txt[:-1]
    mons = []
    for term in txt.split(','):
        if not term:
            continue
        exps = {}
        term = re.sub(r'^[+-]?\d+\*', '', term)
        if re.fullmatch(r'[+-]?\d+', term):
            if int(term) != 0:
                return [()]   # constant -> unit ideal
            continue
        for fac in term.split('*'):
            mm = re.fullmatch(r'u(\d+)(?:\^(\d+))?', fac)
            if not mm:
                raise ValueError('cannot parse msolve term: ' + term[:80])
            i = int(mm.group(1)); k = int(mm.group(2) or 1)
            exps[i] = exps.get(i, 0) + k
        mons.append(tuple(sorted(exps.items())))
    return mons


def hilbert_from_leading(mons, nvars, tmax):
    """h(t) = # monomials of degree t in nvars variables not divisible by any leading monomial, t <= tmax.
    Counts the degree-t multiples of the leading monomials (sparse exponent tuples) with a set."""
    if any(m == () for m in mons):
        return [1] + [0] * tmax, 0
    lead = [m for m in mons]
    degs = [sum(k for _, k in m) for m in lead]
    maxdeg = max(degs) if degs else 0
    h = []
    for t in range(tmax + 1):
        multiples = set()
        for m, dm in zip(lead, degs):
            if dm > t:
                continue
            r = t - dm
            base = dict(m)
            # all monomials of degree r in nvars variables, multiplied onto base
            for c in itertools.combinations_with_replacement(range(1, nvars + 1), r):
                e = dict(base)
                for v in c:
                    e[v] = e.get(v, 0) + 1
                multiples.add(tuple(sorted(e.items())))
        h.append(math.comb(t + nvars - 1, nvars - 1) - len(multiples))
    return h, maxdeg


def hilbert_from_series(q, nvars, tmax):
    """Hilbert function values from the first Hilbert series numerator q (H(t) = sum q_i t^i / (1-t)^nvars)."""
    h = []
    for t in range(tmax + 1):
        s = 0
        for i, qi in enumerate(q):
            if i <= t:
                s += qi * math.comb(t - i + nvars - 1, nvars - 1)
        h.append(s)
    return h


def run_msolve(ms_path, out_path, threads):
    t0 = time.time()
    cmd = ['msolve', '-g', '1', '-t', str(threads), '-f', ms_path, '-o', out_path]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=10 * 3600)
    return {'seconds': round(time.time() - t0, 1), 'returncode': r.returncode, 'stderr_tail': r.stderr[-500:]}


def run_singular(sing_path, out_path, nvars=1):
    t0 = time.time()
    with open(out_path, 'w') as f:
        r = subprocess.run(['Singular', '-q', sing_path], stdout=f, stderr=subprocess.STDOUT, timeout=10 * 3600)
    txt = open(out_path).read()
    res = {'seconds': round(time.time() - t0, 1), 'returncode': r.returncode}
    for key in ['std_time', 'dim', 'vdim', 'hilb1']:
        m = re.search(key + r' = (.*)', txt)
        if m:
            res[key] = m.group(1).strip()
    if 'hilb1' in res:
        try:
            q = [int(x) for x in res['hilb1'].split(',')]
            res['hilbert_function'] = hilbert_from_series(q, nvars, 8)
        except Exception as ex:
            res['hilb_parse_error'] = repr(ex)
    m = re.search(r'dimension \(proj\.\)\s*=\s*(-?\d+).*?degree \(proj\.\)\s*=\s*(\d+)', txt, re.S)
    if m:
        res['proj_dim'], res['proj_degree'] = int(m.group(1)), int(m.group(2))
    res['out_tail'] = txt[-400:]
    return res


def main():
    d = int(sys.argv[1]); parsarg = sys.argv[2]
    pars = {'e': [0], 'o': [1], 'm': [0, 1]}[parsarg]
    pidx = int(sys.argv[sys.argv.index('--pidx') + 1]) if '--pidx' in sys.argv else 0
    threads = int(sys.argv[sys.argv.index('--threads') + 1]) if '--threads' in sys.argv else 16
    tag = sys.argv[sys.argv.index('--tag') + 1] if '--tag' in sys.argv else ''
    p = W.PRIMES30[pidx]
    slice_euler = '--no-slice' not in sys.argv
    job = os.environ.get('W2JOB', f'excl_relax_d{d}{parsarg}')
    base = os.path.join(W.RES, f'excl_relax_d{d}{parsarg}_p{pidx}{tag}')
    F, A, q, rec = W.load()
    Fp, Ap = W.reduce_mod(F, p), W.reduce_mod(A, p)
    t0 = time.time()
    kern = []
    kinfo = []
    for par in pars:
        info, basis = W.e1_kernel(F, A, d, par, p, want_basis=True)
        kinfo.append(info)
        assert info['hilbert_ok'], info
        for Da in (basis or []):
            Na, minfo = W.membership_mod_p(Fp, W.mul(Ap, Da, p), p)
            assert Na is not None, 'kernel element failed membership (should not happen)'
            kern.append((Da, Na))
    print(f'kernel elements: {len(kern)} ({time.time()-t0:.0f}s)', flush=True)
    sysd = W.build_relaxation(F, A, d, kern, pars, p, slice_euler=slice_euler)
    W.save_system(sysd, base + '.system.json')
    W.write_msolve(sysd, base + '.ms', p)
    W.write_singular(sysd, base + '.sing', p)
    nu = len(sysd['unknowns'])
    rec = {'exp': 'relaxation', 'd': d, 'pars': parsarg, 'prime': p, 'unknowns': nu, 'equations': len(sysd['eqs']),
           'n_kernel': len(kern), 'slice_euler': slice_euler, 'parts': sysd['info']['parts'], 'job': job,
           'kernel_info': [{k: v for k, v in ki.items() if k in ('parD', 'kernel mod J', 'hilbert_ok', 'rank J_{d+6}')} for ki in kinfo]}
    if '--no-msolve' not in sys.argv:
        ms = run_msolve(base + '.ms', base + '.ms.out', threads)
        try:
            mons = parse_msolve_leading(base + '.ms.out')
            maxdeg = max((sum(k for _, k in m) for m in mons), default=0)
            # choose tmax: go up to maxdeg + 3 as long as the enumeration cost stays moderate
            degs = [sum(k for _, k in m) for m in mons]
            tmax = 0
            for t in range(1, maxdeg + 4):
                cost = sum(math.comb(t - dm + nu - 1, nu - 1) for dm in degs if dm <= t)
                if cost > 3e7:
                    break
                tmax = t
            h = None
            if tmax >= 1:
                h, _ = hilbert_from_leading(mons, nu, tmax)
            ms.update({'n_leading': len(mons), 'max_leading_degree': maxdeg, 'hilbert_function': h})
            if h is not None:
                ms['proj_empty'] = (0 in h[1:])
        except Exception as ex:
            ms['parse_error'] = repr(ex)
        rec['msolve'] = ms
        print('msolve:', json.dumps(ms)[:600], flush=True)
        W.wlog(exp='relaxation-msolve-partial', d=d, pars=parsarg, prime=p, unknowns=nu, equations=len(sysd['eqs']),
               result=ms.get('proj_empty'), seconds=ms['seconds'], job=job, note='partial record (msolve only)')
    if '--no-singular' not in sys.argv:
        sg = run_singular(base + '.sing', base + '.sing.out', nu)
        rec['singular'] = sg
        print('singular:', json.dumps(sg)[:600], flush=True)
    verdict = 'unknown'
    if rec.get('msolve', {}).get('proj_empty') is True or str(rec.get('singular', {}).get('dim', '')) == '0':
        verdict = 'proj-empty'
    elif 'singular' in rec and 'dim' in rec['singular']:
        verdict = 'affine-dim-' + str(rec['singular']['dim'])
    elif rec.get('msolve', {}).get('proj_empty') is False:
        verdict = 'not-empty-or-undetermined'
    rec['result'] = verdict
    rec['seconds'] = round(time.time() - t0, 1)
    rec['note'] = ('proj-empty => no denominator of degree d with this parity structure over Q-bar other than '
                   'multiples of F (Hilbert-function inequality, prime good in degree 2d)' if verdict == 'proj-empty' else '')
    W.wlog(**rec)


if __name__ == '__main__':
    main()
