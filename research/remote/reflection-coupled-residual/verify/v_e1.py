"""Independent E1-kernel computation by Klein character.
K_d^chi = {D in S_d^chi : A*D in J_{d+6}},   kernel mod J = dim K - dim J_d^chi.
dim K = dim S_d^chi - (rank[J ; A*S_d] - rank J)  in S_{d+6}^{chi*chi_A}.
Also records dim (S/J) in both degrees and compares with the Fermat (all exponents<=4) counts: equality certifies that
rank_p J = rank_Q J, which is the hypothesis under which  ker_Q <= ker_p  (so zeros are rigorous over Q).
Usage: python3 v_e1.py <p> <dmin> <dmax> [chars: e.g. 00,01,10,11]"""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/verify')
from vlib import *
from flint import nmod_mat
p = int(sys.argv[1]); dmin = int(sys.argv[2]); dmax = int(sys.argv[3])
chars = [(int(s[0]), int(s[1])) for s in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['00','01','10','11'])]
Fs, As, qv = load_rec(); F = modp(parse(Fs), p); A = modp(parse(As), p); cA = char_of(A)
Fi = [deriv(F, i, p) for i in range(6)]
cFi = [char_of(f) for f in Fi]
OUT = f'/tmp/hodge/results/reflection-coupled-residual/verify/e1_p{p}.jsonl'
def rank_rows(rows, mons):
    idx = {m: j for j, m in enumerate(mons)}
    return nmod_from_rows(rows, idx, len(mons), p).rank()
for d in range(dmin, dmax + 1):
    for c in chars:
        t0 = time.time()
        cT = chi_mul(c, cA)
        monsD = monos_chi(d, c); monsT = monos_chi(d + 6, cT)
        # J rows in degree d+6 with character cT: m*F_i, m of degree d+1 with character cT*cFi
        Jrows = []
        for i in range(6):
            for m in monos_chi(d + 1, chi_mul(cT, cFi[i])):
                Jrows.append(shift(Fi[i], m))
        Arows = [shift(A, m) for m in monsD]
        rJ = rank_rows(Jrows, monsT)
        rJA = rank_rows(Jrows + Arows, monsT)
        # J in degree d with character c
        Jd = []
        for i in range(6):
            if d - 5 >= 0:
                for m in monos_chi(d - 5, chi_mul(c, cFi[i])):
                    Jd.append(shift(Fi[i], m))
        rJd = rank_rows(Jd, monsD) if Jd else 0
        dimK = len(monsD) - (rJA - rJ)
        ker_modJ = dimK - rJd
        quotD = len(monsD) - rJd; quotT = len(monsT) - rJ
        fD = fermat_count(d, c); fT = fermat_count(d + 6, cT)
        rec = {'p': p, 'd': d, 'chi_D': list(c), 'chi_T': list(cT), 'dim_S_d_chi': len(monsD), 'dim_S_T': len(monsT),
               'rank_J_T': rJ, 'rank_JA': rJA, 'rank_J_d': rJd, 'dim_RD': quotD, 'dim_RT': quotT,
               'fermat_RD': fD, 'fermat_RT': fT, 'hilbert_ok': (quotD == fD and quotT == fT),
               'dimK': dimK, 'ker_modJ': ker_modJ, 'secs': round(time.time() - t0, 1)}
        print(json.dumps(rec), flush=True)
        with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
