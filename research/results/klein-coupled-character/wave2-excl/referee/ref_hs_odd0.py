"""Hilbert-function sandwich for (d=7, Y-odd, char 0) from a msolve -g 1 leading ideal file (24 variables):
HF(S/I(L)) <= HF(S/I_Q) <= HF(S/I_p) = HF(S/in(I_p)); L = cone over Segre(C x P^2), C = CI of three quadrics in P^4,
so HF_L(e) = HF_C(e) * C(e+2,2) with HF_C(e) = coefficient of t^e in (1+t)^3/(1-t)^2 (= 8e-4 for e >= 2).
usage: python3 ref_hs_odd0.py <leadfile> <tag>"""
import sys, re, subprocess, time
from math import comb
leadf, tag = sys.argv[1], sys.argv[2]
txt = open(leadf).read(); body = txt[txt.index('[') + 1: txt.rindex(']')]
items = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
nv = 24
sg = f'ring R = 32003, ({",".join(f"u{i}" for i in range(nv))}), dp;\nideal LEAD = ' + ',\n'.join(items) + ';\nideal SL = std(LEAD);\n"n_lead", size(LEAD), "gb", size(SL), "dim", dim(SL), "mult", mult(SL);\nintvec h = hilb(SL, 1);\n"HS1", h;\nquit;\n'
open(f'hs_{tag}.sing', 'w').write(sg)
t = time.time(); r = subprocess.run(['Singular', '-q', f'hs_{tag}.sing'], capture_output=True, text=True, timeout=7200)
out = r.stdout; print(out[-600:], r.stderr[:200], f'[{time.time()-t:.0f}s]', flush=True)
m = re.search(r'HS1\s+([\d\-, \n]+)', out); num = [int(x) for x in re.split(r'[,\s]+', m.group(1).strip()) if x not in ('', '-')]
while num and num[-1] == 0: num.pop()
def hf(num, n, e): return sum(num[i] * comb(e - i + n - 1, n - 1) for i in range(min(e, len(num) - 1) + 1))
HFP = [hf(num, nv, e) for e in range(0, 16)]
def hfC(e): return 1 if e == 0 else (5 if e == 1 else 8*e - 4)
HFL = [hfC(e) * comb(e + 2, 2) for e in range(0, 16)]
print('HF(S/in I_p) e=0..15:', HFP); print('HF_L (Segre)  e=0..15:', HFL)
agree_from = next((e for e in range(16) if all(HFP[f] == HFL[f] for f in range(e, 16))), None)
print('agree for all e >=', agree_from, '(sandwich closes iff agreement persists: both are polynomials of degree 3 for e >= 2, so agreement on 4 consecutive values e>=2 suffices)')
print('HS_ODD0_DONE', agree_from is not None and agree_from <= 12)
