"""d=8 eigen sectors: msolve GB + radical membership against the lifted family S_2(sector-compatible quadrics).C.
usage: python3 w2_radical8.py <par> <k> <threads>"""
import sys, time, json, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
par, k, thr = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
F, A, q = load(); p = 1000003; ctx = Ctx(F, A, p)
OUT = 'local_radical8'; os.makedirs(OUT, exist_ok=True)
tag = f'd8_{"even" if par == 0 else "odd"}_k{k}'
t0 = time.time()
s = hess_system(ctx, 8, [par], [k])
ms = f'{OUT}/{tag}.ms'; write_msolve(s, ms)
lead = f'{OUT}/{tag}.lead'
r = subprocess.run(['msolve', '-f', ms, '-o', lead, '-g', '1', '-t', str(thr)], capture_output=True, text=True, timeout=36000)
mons = parse_msolve_leading(lead); irr, info = leading_ideal_irrelevant(mons, s['n'])
print(f'{tag}: n={s["n"]} eqs={s["n_eqs"]} leading ideal size {len(mons)} irrelevant={irr} [{time.time()-t0:.0f}s]', flush=True)
if irr:
    log({'exp': 'hessian-relaxation-eigen', 'family': tag, 'd': 8, 'par': par, 'k': k, 'prime': p, 'result': 'EMPTY (leading ideal irrelevant: no sigma-eigen denominator of this parity/character at degree 8)', 'seconds': round(time.time()-t0, 1), 'job': 'local', 'note': f'n_vars {s["n"]} n_eqs {s["n_eqs"]}'})
    print('RADICAL_DONE empty'); sys.exit(0)
gbf = f'{OUT}/{tag}.gb'
r = subprocess.run(['msolve', '-f', ms, '-o', gbf, '-g', '2', '-t', str(thr)], capture_output=True, text=True, timeout=36000)
gb = msolve_gb_to_singular(gbf); ngb = gb.count(',') + 1
Cq = [{tuple(int(x) for x in kk.split(',')): Fraction(v) for kk, v in qd.items()} for qd in json.load(open('local_d6/C_quadrics_exact.json'))]
# multiplier quadrics: parity par (D parity) since D_lambda is even; character k (D_lambda char 0)
mult = wmonomials(2, par, k)
lin, minors, cq, Mindex = lifted_family_generators_general(s['vars'], mult, Cq)
print(f'  multipliers {len(mult)}: {mult}; lifted family: {len(Mindex)} matrix coords, {len(lin)} vanishing, {len(minors)} minors, {len(cq)} C-quadric images; GB size {ngb}', flush=True)
sg = singular_radical_membership(None, s, {'vanishing-coordinates': lin, 'rank-one-minors': minors, 'C-quadrics': cq}, f'{OUT}/{tag}.sing', p, maxk=8, from_msolve_gb=gb)
t = time.time()
r = subprocess.run(['Singular', '-q', sg], capture_output=True, text=True, timeout=36000)
out = r.stdout.strip(); ok = out.endswith('total failures 0')
print(out[-1200:], r.stderr[:300], f'[{time.time()-t:.0f}s]', flush=True)
log({'exp': 'radical-membership-lifted', 'family': tag, 'd': 8, 'par': par, 'k': k, 'prime': p,
     'result': 'V(relaxation) subset of lifted family S_2.C' if ok else 'FAILURE: extra components beyond S_2.C', 'seconds': round(time.time()-t0, 1), 'job': 'local',
     'note': f'n_vars {s["n"]}, n_eqs {s["n_eqs"]}, GB size {ngb}, multipliers {len(mult)}; ' + out[-300:].replace('\n', ' | ')})
print('RADICAL_DONE', ok)
