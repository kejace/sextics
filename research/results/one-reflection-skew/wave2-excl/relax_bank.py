"""Build the bank of relaxed correction columns B_{D_a}(K) for all family directions a and all Koszul K (all characters),
grouped by row sector psi of S_{2d}; plus the F*S_{2d-6} columns.  Saved as d<d>_bank.pkl:
{psi: {'rows': nrows, 'F': [sparse cols], 'corr': [(a, sparse col)]}}.  Usage: python3 relax_bank.py <d>"""
import sys, pickle, time
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx
d = int(sys.argv[1]); n2 = 2 * d
obj = pickle.load(open(f'd{d}_fam.pkl', 'rb'))
items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
fam = w2x.Family(d, items)
t0 = time.time()
bank = {psi: {'rows': len(w2x.secmons(n2, psi)), 'F': [], 'corr': []} for psi in w2x.ALLCHR}
for m in w2x.mons(n2 - 6):
    P = w2x.F7 * w2x.monom(m); psi = w2x.chr_add(w2x.chr_of(m), w2x.CHI_F)
    bank[psi]['F'].append(w2x.vec_of(P, n2, psi))
nc = 0
for a in range(fam.n):
    _, Da, _, _, chia, _ = fam.items[a]
    for psi in w2x.ALLCHR:
        chiK = w2x.chr_add(psi, w2x.chr_neg(chia), w2x.chr_neg(w2x.CHI_A))
        for (i, j, m) in w2x.koszul_cols(d, chiK):
            P = w2x.BD_koszul(Da, i, j, m)
            if P.is_zero(): continue
            bank[psi]['corr'].append((a, w2x.vec_of(P, n2, psi))); nc += 1
    if a % 10 == 0: print('a', a, 'cols', nc, round(time.time() - t0, 1), flush=True)
pickle.dump(bank, open(f'd{d}_bank.pkl', 'wb'))
print('bank done', nc, 'correction columns', round(time.time() - t0, 1), flush=True)
