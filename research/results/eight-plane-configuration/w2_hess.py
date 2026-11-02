"""Wave-2 job: exact (mod p) Hessian-relaxation quadrics for the Jacobian family.

For D = w(F) (w a vector of forms of degree d-5) the E2 identity reduced modulo J(F) is
    A * Hess_F(w, w) in J(F)_{2d}   <=>   [Hess_F(w,w)] in Ann_R(A)_{2d-6}.
This is a system of quadrics in the w-coefficients.  For each (d, parity) requested we compute the
quadrics mod p, their rank, the Macaulay ranks in low degree (decides emptiness of the projective
variety when the rank reaches the number of monomials), write an msolve input file, and optionally
run msolve (grevlex Groebner basis) with a time limit.

usage: python3 w2_hess.py --degrees 6,7 --pars 0,1 [--macaulay-max 3] [--msolve 1] [--msolve-timeout 7200]
                          [--threads 8] [--tag hess]
"""
import sys, os, json, time, argparse, subprocess, itertools
ap = argparse.ArgumentParser()
ap.add_argument('--degrees', default='6')
ap.add_argument('--pars', default='0,1')
ap.add_argument('--macaulay-max', type=int, default=3)
ap.add_argument('--macaulay-maxcols', type=int, default=60000)
ap.add_argument('--msolve', type=int, default=0)
ap.add_argument('--msolve-timeout', type=int, default=7200)
ap.add_argument('--threads', type=int, default=8)
ap.add_argument('--tag', default='hess')
args = ap.parse_args()
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
from flint import nmod_mat

OUT = f'{WD}/{args.tag}.jsonl'


def nmon(n, k):
    from math import comb
    return comb(n + k - 1, k)


def msolve_leading_check(gb_text, nW):
    """From msolve -g 2 output (list of polynomials), check whether pure powers of all variables occur as
    leading monomials (grevlex, msolve prints terms in decreasing order): if so the projective variety is empty."""
    import re
    polys = [s.strip() for s in gb_text.replace('\n', '').split(',') if s.strip()]
    pure = set()
    nlead = 0
    for s in polys:
        s = s.strip('[] ')
        if not s:
            continue
        lead = s.split('+')[0].split('-')[0] if not s.startswith('-') else s[1:].split('+')[0].split('-')[0]
        lead = lead.strip()
        nlead += 1
        m = re.fullmatch(r'(?:\d+\*)?w(\d+)(?:\^(\d+))?', lead)
        if m:
            pure.add(int(m.group(1)))
    return {'n_gb': nlead, 'pure_power_vars': len(pure), 'all_pure': len(pure) == nW}


for d in [int(x) for x in args.degrees.split(',')]:
    for pD in [int(x) for x in args.pars.split(',')]:
        t0 = time.time()
        res = hess_quadrics(d, pD, p=PRIME, verbose=True)
        rows = quadric_rows(res)
        nW = res['nW']
        rec = {'event': 'hess', 'd': d, 'par': pD, 'p': PRIME, 'nW': nW, 'npairs': res['npairs'],
               'dimR_2d_odd': res['dimR_2d_odd'], 'rank_quadrics': res['rank_quadrics'],
               'expected_conditions_note': 'rank <= min(dim Sym^2 W, rank of A: R_{2d-6}^+ -> R_{2d}^-)',
               't_quadrics': round(time.time() - t0, 1)}
        # Macaulay ranks: degree k matrix has rows = rank_quadrics * nmon(nW, k-2), cols = nmon(nW, k)
        mac = {}
        for k in range(2, args.macaulay_max + 1):
            cols = nmon(nW, k)
            if cols > args.macaulay_maxcols:
                mac[str(k)] = {'skipped': True, 'cols': cols}
                break
            t1 = time.time()
            rk, R, C = macaulay_rank(rows, nW, k, PRIME)
            mac[str(k)] = {'rank': rk, 'rows': R, 'cols': C, 'full': rk == C, 'codim_hilbert': C - rk,
                           'seconds': round(time.time() - t1, 1)}
            print(f'[macaulay] d={d} par={pD} deg={k} rank={rk} of {C} ({R} rows) -> Hilbert function {C-rk}',
                  file=sys.stderr, flush=True)
            if rk == C:
                rec['projective_variety_empty_modp'] = True
                break
            if rk == C - 1 and k >= 2:
                rec['hilbert_1_at_degree'] = k   # a single (reduced) point
        rec['macaulay'] = mac
        path = f'{WD}/{args.tag}-quadrics-d{d}-p{pD}.ms'
        nrows = write_msolve(res, path)
        rec['msolve_file'] = path
        # save the quadrics (reduced rows) for later use
        with open(f'{WD}/{args.tag}-quadrics-d{d}-p{pD}.json', 'w') as f:
            json.dump({'d': d, 'par': pD, 'p': PRIME, 'nW': nW, 'Wm': [[i, list(m)] for i, m in res['Wm']],
                       'rows': [[[list(k), v] for k, v in r.items()] for r in rows]}, f)
        if args.msolve and not rec.get('projective_variety_empty_modp'):
            t1 = time.time()
            outp = f'{WD}/{args.tag}-gb-d{d}-p{pD}.txt'
            cmd = ['msolve', '-g', '2', '-t', str(args.threads), '-v', '1', '-f', path, '-o', outp]
            try:
                pr = subprocess.run(cmd, capture_output=True, text=True, timeout=args.msolve_timeout)
                rec['msolve'] = {'returncode': pr.returncode, 'seconds': round(time.time() - t1, 1),
                                 'stderr_tail': pr.stderr[-800:], 'stdout_tail': pr.stdout[-400:]}
                if os.path.exists(outp):
                    txt = open(outp).read()
                    rec['msolve'].update(msolve_leading_check(txt, nW))
                    rec['msolve']['gb_file'] = outp; rec['msolve']['gb_bytes'] = len(txt)
            except subprocess.TimeoutExpired:
                rec['msolve'] = {'timeout': args.msolve_timeout}
        rec['seconds'] = round(time.time() - t0, 1)
        jlog(OUT, rec)
        print(json.dumps({k: v for k, v in rec.items() if k != 'macaulay'}), flush=True)
        print(json.dumps(rec['macaulay']), flush=True)
