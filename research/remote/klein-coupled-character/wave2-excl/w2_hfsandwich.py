"""Hilbert-function sandwich: rigorous V(I_Q) = L from a mod-p Groebner basis.
  HF(S/I(L)) <= HF(S/rad I_Q) <= HF(S/I_Q) <= HF(S/I_p)   (I(L) exact ideal of the lifted family, I_p the mod-p relaxation ideal)
If the Hilbert series of S/in(I_p) (from the msolve leading ideal) and of S/I(L) (exact, Singular) differ by a polynomial
(i.e. HF agree in all large degrees) then rad(I_Q) = I(L) in large degrees, hence V(I_Q) = L over Qbar.
usage: python3 w2_hfsandwich.py <tag> <d> <par> <chars|all> <leading-ideal-file (msolve -g 1 or -g 2 output)>"""
import sys, time, json, subprocess, os, re
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import *
tag, d, par, chars, leadf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
chars = list(range(7)) if chars == 'all' else [int(c) for c in chars.split(',')]
F, A, q = load(); p = 1000003; ctx = Ctx(F, A, p)
OUT = 'local_hfsandwich'; os.makedirs(OUT, exist_ok=True)
s = hess_system(ctx, d, [par], chars)
n = s['n']; names = ','.join(f'u{a}' for a in range(n))
# leading monomials (from -g 1 file, or from the leading terms of a -g 2 GB file)
txt = open(leadf).read()
body = txt[txt.index('[') + 1: txt.rindex(']')]
items = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
if '-g 2' in txt or any('+' in it for it in items):
    # full GB: leading term = first term (msolve prints terms in decreasing grevlex order)
    lead = []
    for it in items:
        first = re.split(r'(?<!\^)[+-](?=\d)', it)[0] if it[0] not in '+-' else re.split(r'(?<!\^)[+-](?=\d)', it[1:])[0]
        m = '*'.join(f for f in first.split('*') if not re.fullmatch(r'\d+', f))
        lead.append(m if m else '1')
else:
    lead = items
lead = [m for m in lead]
Cq = [{tuple(int(x) for x in kk.split(',')): Fraction(v) for kk, v in qd.items()} for qd in json.load(open('local_d6/C_quadrics_exact.json'))]
mult = wmonomials(d - 6, par, None)
mult = [m for m in mult if (wt(m) % 7) in chars] if len(chars) < 7 else mult
lin, minors, cq, Mindex = lifted_family_generators_general(s['vars'], mult, Cq)
# exact ideal of L via elimination: variables M (matrix coords) + lambda_1..5 + mu_j (one per multiplier)
lam = [f'l{i}' for i in range(1, 6)]; mu = [f'm{j}' for j in range(len(mult))]
elim_vars = ','.join(lam + mu)
param = []
for (i, j), a in Mindex.items():
    param.append(f'u{a}-l{i}*m{j}')
import math
cqs = []
for qd in Cq:
    den = 1
    for v in qd.values(): den = den * v.denominator // math.gcd(den, v.denominator)
    cqs.append('+'.join(f'({int(v*den)})*l{i}*l{ip}' for (i, ip), v in sorted(qd.items())))
sg = f'''LIB "poly.lib";
ring R = 0, ({names},{elim_vars}), dp;
ideal P = {','.join(param + cqs)};
ideal IL = eliminate(P, {'*'.join(lam + mu)});
ring S = 0, ({names}), dp;
ideal IL = imap(R, IL);
ideal ILfull = IL, {','.join(lin) if lin else '0'};
ideal SIL = std(ILfull);
"IL_gb_size", size(SIL), "dim", dim(SIL), "mult", mult(SIL);
intvec h1 = hilb(SIL, 1);
"HS_L", h1;
ideal LEAD = {','.join(lead)};
ideal SLEAD = std(LEAD);
"lead_dim", dim(SLEAD), "lead_mult", mult(SLEAD);
intvec h2 = hilb(SLEAD, 1);
"HS_P", h2;
quit;
'''
fn = f'{OUT}/{tag}.sing'; open(fn, 'w').write(sg)
t = time.time()
r = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=7200)
out = r.stdout
print(out[-1500:], r.stderr[:300], f'[{time.time()-t:.0f}s]', flush=True)
def vec(name):
    m = re.search(name + r'\s+([\d\-, ]+)', out)
    return [int(x) for x in re.split(r'[,\s]+', m.group(1).strip()) if x not in ('', '-')] if m else None
hL, hP = vec('HS_L'), vec('HS_P')
# Singular's hilb(..,1) intvec: coefficients of the first Hilbert series numerator (last entry 0 marker)
def strip(v): 
    v = v[:]
    while v and v[-1] == 0: v.pop()
    return v
hL, hP = strip(hL), strip(hP)
diff = [a - b for a, b in zip(hL + [0] * (len(hP)), hP + [0] * (len(hL)))]
# check (1-t)^n divides diff: repeated synthetic division by (1-t) i.e. polynomial in t: divide by (1 - t) n times
def div_by_1_minus_t(c):
    # c coefficients of polynomial; divide by (1 - t): q with (1-t) q = c  => q_k = c_k + q_{k-1}
    qc = []; acc = 0
    for x in c:
        acc = x + acc; qc.append(acc)
    rem = qc[-1] if qc else 0
    return qc[:-1], rem
cur = diff[:]; divisible = True
for _ in range(n):
    cur, rem = div_by_1_minus_t(cur)
    if rem != 0: divisible = False; break
verdict = 'RIGOROUS: HF(S/I_p) and HF(S/I(L)) agree in all large degrees => V(I_Q) = L over Qbar' if divisible else 'NOT CONCLUSIVE: Hilbert functions differ in large degrees (I_p not reduced enough, or extra components)'
print('HS_L', hL, 'HS_P', hP, 'divisible by (1-t)^n:', divisible, '->', verdict, flush=True)
log({'exp': 'hilbert-sandwich', 'family': tag, 'd': d, 'par': par, 'k': chars if len(chars) < 7 else 'all', 'prime': p, 'result': verdict, 'job': 'local',
     'note': f'n={n}; first Hilbert series numerators: L {hL}, in(I_p) {hP}; dims/mults in local_hfsandwich/{tag}.out'})
open(f'{OUT}/{tag}.out', 'w').write(out)
print('HF_DONE', divisible)
