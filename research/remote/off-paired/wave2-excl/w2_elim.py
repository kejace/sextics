"""EXACT sector exclusion via Singular over Q: elimination of the Koszul coordinates.

Sector (d, c).  Rational data: D_a (basis of the colon sector modulo the Euler gauge; verified N_a(F) = A D_a exactly),
Koszul kernel S_k (verified S_k(F) = 0 exactly), exact B products (fmpq_mpoly), exact normal forms modulo F*S_{2d-6}[sector]
(rational rref).  The bilinear ideal I in Q[t, s] (lossless for Euler-normalized certificates with D in the sector) is
eliminated over Q:  E = I  ∩  Q[t]  (Singular `eliminate`, exact).  E is homogeneous; `dim(std(E)) == 0` (affine cone = {0})
means every certificate has t = 0, i.e. D in the gauge F*S_{d-6} (illegal).  No modular step is involved in the conclusion.
usage: python3 w2_elim.py --d 7 --chars 00,01 --job NAME [--timeout S] [--prime P (only to find the bases)]
"""
import sys, os, json, time, argparse, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib, w2_bilinear as WB, w2_cert as WC
from w2lib import char_of, cadd, monos_char
from hodge.core import NV, shift, reduce_mod
from fractions import Fraction
from flint import fmpq_mat, fmpq

class RatReducer:
    """Exact normal-form coordinates modulo the Q-row space of polys within mons (rref over Q)."""
    def __init__(self, polys, mons):
        self.mons = mons; self.idx = {m: i for i, m in enumerate(mons)}
        R, C = len(polys), len(mons)
        if R == 0:
            self.pivots = []; self.rows = []; self.std = list(range(C)); return
        M = fmpq_mat(R, C)
        for r, P in enumerate(polys):
            for e, v in P.items(): M[r, self.idx[e]] = fmpq(v.numerator, v.denominator)
        Mr, rank = M.rref()
        piv = []; rows = []
        for i in range(rank):
            j = 0
            while j < C and Mr[i, j] == 0: j += 1
            piv.append(j)
            rows.append({k: Fraction(int(Mr[i, k].p), int(Mr[i, k].q)) for k in range(j, C) if Mr[i, k] != 0})
        self.pivots = piv; self.rows = rows
        ps = set(piv); self.std = [j for j in range(C) if j not in ps]
    def coords(self, P):
        v = {self.idx[e]: Fraction(x) for e, x in P.items()}
        for j, row in zip(self.pivots, self.rows):
            c = v.get(j, 0)
            if c:
                for k, x in row.items():
                    nv = v.get(k, 0) - c * x
                    if nv: v[k] = nv
                    else: v.pop(k, None)
        sidx = {j: i for i, j in enumerate(self.std)}
        return {sidx[j]: x for j, x in v.items() if x}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True); ap.add_argument('--chars', required=True); ap.add_argument('--job', required=True)
    ap.add_argument('--timeout', type=int, default=1800); ap.add_argument('--prime', type=int, default=WB.PM); ap.add_argument('--tsub', default='all')
    a = ap.parse_args()
    F, A, q = w2lib.load(); cA = char_of(A); p = a.prime
    chars = [(x, y) for x in range(6) for y in range(6)] if a.chars == 'all' else [(int(s[0]), int(s[1])) for s in a.chars.split(',')]
    outdir = f'{w2lib.PK}/elim_{a.job}'; os.makedirs(outdir, exist_ok=True)
    recfile = f'{w2lib.PK}/receipts_{a.job}.jsonl'
    for c in chars:
        t0 = time.time()
        sysd, info = WB.build_fast(F, A, q, a.d, c, p, 'sector', a.tsub, addJ=False, verbose=False)
        base = {'exp': 'elimQ', 'job': a.job, 'd': a.d, 'char': c, 'tsub': a.tsub, 'family': f'all D in sector ({a.tsub})', 'nt': info.get('nt'), 'ns': info.get('ns')}
        def emit(rec):
            rec['seconds'] = round(time.time() - t0, 1)
            w2lib.log(rec); open(recfile, 'a').write(json.dumps(rec, default=str) + '\n')
            print(json.dumps({k: rec.get(k) for k in ('d', 'char', 'nt', 'ns', 'neq', 'elim_gens', 'dim', 'result', 'seconds')}, default=str), flush=True)
        if sysd is None:
            emit(dict(base, result='PROVEN-empty-sector(colon dim 0 mod p bounds the rational colon dim)')); continue
        try:
            Dq, Nq, Sq = WC.rational_data(F, A, q, sysd['Ds'], sysd['Ns'], sysd['syz'], p)
        except ValueError as ex:
            emit(dict(base, result=f'rational-data-failed: {ex}')); continue
        nt, ns = len(Dq), len(Sq)
        EB = WC.ExactB(Dq, Nq, Sq)
        sector = cadd(cA, cadd(c, c)); mons2d = monos_char(2 * a.d, sector); mons_low = monos_char(2 * a.d - 6, sector)
        red = RatReducer([shift(F, m, 1) for m in mons_low], mons2d)
        eqs = {}
        def addeq(coords, key):
            for r, v in coords.items():
                eqs.setdefault(r, {}); eqs[r][key] = eqs[r].get(key, 0) + v
        for aa in range(nt):
            for bb in range(nt):
                addeq(red.coords(EB.B(aa, 't', bb)), (min(aa, bb), max(aa, bb)))
            for k in range(ns):
                addeq(red.coords(EB.B(aa, 's', k)), (aa, nt + k))
        eqs = [e for e in eqs.values() if any(e.values())]
        names = [f't{i}' for i in range(nt)] + [f's{k}' for k in range(ns)]
        def vn(i): return f't{i}' if i < nt else f's{i - nt}'
        polys = []
        for e in eqs:
            terms = []
            for (i, j), v in e.items():
                if v == 0: continue
                cs = f'({v.numerator}/{v.denominator})' if v.denominator != 1 else f'({v.numerator})'
                terms.append(f'{cs}*{vn(i)}*{vn(j)}')
            polys.append('+'.join(terms))
        sing = f'''ring r = 0, ({",".join(names)}), dp;
ideal I = {",\n".join(polys)};
poly ps = {"*".join(f"s{k}" for k in range(ns)) if ns else "1"};
ideal E = eliminate(I, ps);
ideal G = std(E);
int dm = dim(G);
print("ELIM_GENS " + string(size(E)));
print("DIM " + string(dm));
print("LEAD " + string(lead(G)));
quit;
'''
        path = f'{outdir}/d{a.d}_c{c[0]}{c[1]}.sing'; open(path, 'w').write(sing)
        try:
            r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=a.timeout)
            out = r.stdout
        except subprocess.TimeoutExpired:
            emit(dict(base, neq=len(eqs), result='singular-timeout')); continue
        m1 = re.search(r'ELIM_GENS (\d+)', out); m2 = re.search(r'DIM (-?\d+)', out)
        if not m2:
            emit(dict(base, neq=len(eqs), result='singular-failed', stderr=(r.stderr or out)[-400:])); continue
        dm = int(m2.group(1))
        rec = dict(base, neq=len(eqs), elim_gens=int(m1.group(1)) if m1 else None, dim=dm, lead=out.split('LEAD ')[-1][:300] if 'LEAD ' in out else None)
        # dim is taken in the ring with the s variables, which E does not involve: dim_t = dm - ns
        rec['dim_t'] = dm - ns if dm >= 0 else -1
        if dm - ns <= 0:
            rec['result'] = 'PROVEN-no-denominator-in-sector(exact elimination over Q: E is m-primary)'
        else:
            rec['result'] = f'elimination ideal has positive dimension {dm - ns} (not excluded by this method)'
        emit(rec)

if __name__ == '__main__':
    main()
