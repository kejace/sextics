"""Radical-membership test V(I_relax) subset L (lifted family) for a relaxation system: msolve -g 2 GB mod p, then Singular reduce.
usage: python3 w2_radical.py <tag> <d> <par> <chars:comma or all> <threads>"""
import sys, time, json, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
tag, d, par, chars, thr = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
chars = list(range(7)) if chars == 'all' else [int(c) for c in chars.split(',')]
F, A, q = load(); p = 1000003; ctx = Ctx(F, A, p)
OUT = 'local_radical'; os.makedirs(OUT, exist_ok=True)
t0 = time.time()
s = hess_system(ctx, d, [par], chars)
ms = f'{OUT}/{tag}.ms'; write_msolve(s, ms)
gbf = f'{OUT}/{tag}.gb'
if not (os.path.exists(gbf) and os.path.getsize(gbf) > 0):
    r = subprocess.run(['msolve', '-f', ms, '-o', gbf, '-g', '2', '-t', str(thr)], capture_output=True, text=True, timeout=36000)
    print('msolve rc', r.returncode, r.stderr[-200:], f'[{time.time()-t0:.0f}s]', flush=True)
gb = msolve_gb_to_singular(gbf)
ngb = gb.count(',') + 1
Cq = [{tuple(int(x) for x in k.split(',')): Fraction(v) for k, v in qd.items()} for qd in json.load(open('local_d6/C_quadrics_exact.json'))]
lin, minors, cq, Mindex = lifted_family_generators(s['vars'], par, Cq)
print(f'{tag}: n={s["n"]} eqs={s["n_eqs"]} GB size {ngb}; lifted family: {len(Mindex)} matrix coords, {len(lin)} vanishing coords, {len(minors)} minors, {len(cq)} C-quadric images', flush=True)
sg = singular_radical_membership(None, s, {'vanishing-coordinates': lin, 'rank-one-minors': minors, 'C-quadrics': cq}, f'{OUT}/{tag}.sing', p, maxk=8, from_msolve_gb=gb)
t = time.time()
r = subprocess.run(['Singular', '-q', sg], capture_output=True, text=True, timeout=36000)
out = r.stdout.strip()
print(out[-1500:], r.stderr[:300], f'[{time.time()-t:.0f}s]', flush=True)
ok = out.strip().endswith('total failures 0')
log({'exp': 'radical-membership-lifted', 'family': tag, 'd': d, 'par': par, 'k': chars if len(chars) < 7 else 'all', 'prime': p,
     'result': 'V(relaxation) subset of lifted family S_1.C (all generators of I_L in the radical)' if ok else 'FAILURE: some generator of I_L not in radical (extra components)',
     'seconds': round(time.time()-t0, 1), 'job': 'local', 'note': f'n_vars {s["n"]}, n_eqs {s["n_eqs"]}, GB size {ngb}; ' + out[-400:].replace('\n', ' | ')})
print('RADICAL_DONE', ok)
