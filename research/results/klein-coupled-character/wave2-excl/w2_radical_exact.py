"""EXACT (characteristic 0) radical-membership test V(I_relax) subset L for an eigen sector at d = 7 (or 6):
builds the relaxation system over Q (hess_system_exact), computes the reduced GB over Q with msolve (char 0), then
tests g^k in I over Q with Singular for every generator g of the lifted-family ideal I_L.
usage: python3 w2_radical_exact.py <tag> <d> <par> <chars:comma or all> <threads>"""
import sys, time, json, subprocess, os, math
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
tag, d, par, chars, thr = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
chars = list(range(7)) if chars == 'all' else [int(c) for c in chars.split(',')]
F, A, q = load()
OUT = 'local_radical_exact'; os.makedirs(OUT, exist_ok=True)
t0 = time.time()
s = hess_system_exact(F, A, d, [par], chars)
n = s['n']
# msolve input over Q (integer-scaled equations)
names = ','.join(f'u{a}' for a in range(n))
polys = []
for eq in s['eqs']:
    den = 1
    for v in eq.values(): den = den * v.denominator // math.gcd(den, v.denominator)
    polys.append('+'.join((f'{int(c*den)}*u{a}*u{b}' if a != b else f'{int(c*den)}*u{a}^2') for (a, b), c in sorted(eq.items())))
ms = f'{OUT}/{tag}.ms'
polys = [pl.replace('+-', '-') for pl in polys]
open(ms, 'w').write(names + '\n0\n' + ',\n'.join(polys) + '\n')
gbf = f'{OUT}/{tag}.gb'
r = subprocess.run(['msolve', '-f', ms, '-o', gbf, '-g', '2', '-t', str(thr)], capture_output=True, text=True, timeout=36000)
print('msolve (char 0) rc', r.returncode, r.stderr[-300:], f'[{time.time()-t0:.0f}s]', flush=True)
gb = msolve_gb_to_singular(gbf); ngb = gb.count(',') + 1
if ngb <= 1:
    print('UNEXPECTED trivial GB from msolve char 0 - falling back to Singular std over Q', flush=True); gb = None
Cq = [{tuple(int(x) for x in kk.split(',')): Fraction(v) for kk, v in qd.items()} for qd in json.load(open('local_d6/C_quadrics_exact.json'))]
mult = wmonomials(d - 6, par, None)   # multipliers of degree d-6 with the parity of D (D_lambda is even, char 0): char must match the sector
mult = [m for m in mult if (wt(m) % 7) in chars] if len(chars) < 7 else mult
lin, minors, cq, Mindex = lifted_family_generators_general(s['vars'], mult, Cq)
print(f'{tag}: n={n} eqs={s["n_eqs"]} exact GB size {ngb}; multipliers {len(mult)}; lifted family: {len(Mindex)} coords, {len(lin)} vanishing, {len(minors)} minors, {len(cq)} C-quadric images', flush=True)
body = f'ring R = 0, ({names}), dp;\n' + ((f'ideal SI = ' + gb + ';\nattrib(SI, "isSB", 1);\n') if gb else ('ideal I = ' + ',\n'.join(polys) + ';\nideal SI = std(I);\n')) + '"gb_size", size(SI), "dim", dim(SI);\nint k; int bad = 0;\n'
for name, gens in {'vanishing-coordinates': lin, 'rank-one-minors': minors, 'C-quadrics': cq}.items():
    body += f'"group {name}: {len(gens)} generators";\n'
    for g in gens:
        body += f'k = 1; while (k <= 8) {{ if (reduce(({g})^k, SI) == 0) {{ break; }} k = k + 1; }} if (k > 8) {{ bad = bad + 1; "FAIL", "{g[:60]}"; }}\n'
body += '"total failures", bad;\nquit;\n'
sg = f'{OUT}/{tag}.sing'; open(sg, 'w').write(body)
t = time.time()
r = subprocess.run(['Singular', '-q', sg], capture_output=True, text=True, timeout=36000)
out = r.stdout.strip(); ok = out.endswith('total failures 0')
print(out[-1200:], r.stderr[:300], f'[{time.time()-t:.0f}s]', flush=True)
log({'exp': 'radical-membership-lifted-EXACT', 'family': tag, 'd': d, 'par': par, 'k': chars if len(chars) < 7 else 'all', 'prime': 0,
     'result': 'RIGOROUS over Q: V(relaxation) subset of the lifted family S_{d-6}.C (every generator of I_L has a power in the exact relaxation ideal)' if ok else 'FAILURE over Q: extra components beyond the lifted family',
     'seconds': round(time.time()-t0, 1), 'job': 'local', 'note': f'n_vars {n}, n_eqs {s["n_eqs"]}, exact GB size {ngb}; ' + out[-300:].replace('\n', ' | ')})
print('RADICAL_EXACT_DONE', ok)
