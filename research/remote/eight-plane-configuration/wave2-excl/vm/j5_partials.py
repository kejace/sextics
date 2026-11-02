"""(a)+(b) d=5: Koszul-reduced systems for the six partials (two primes + exact fmpz ranks)."""
import sys, time, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl/vm')
from kred import *
P1, P2 = 2147483647, 2147483587
for i in range(6):
    w = unit_field(i)
    S = KoszulSystem(w, parity_rows=None, verbose=True)
    r1 = S.solve_mod(P1); r2 = S.solve_mod(P2)
    rM, tM = S.exact_rank(False); rMb, tMb = S.exact_rank(True)
    rec = {'exp': 'kred-partial', 'D': f'dF/d{core.VARS[i]}', 'd': 5, 'rows': S.nrows, 'cols': S.ncols,
           'p1': P1, 'consistent_p1': r1['consistent'], 'rank_p1': r1['rank_p'], 'rank_aug_p1': r1['rank_aug_p'],
           'p2': P2, 'consistent_p2': r2['consistent'], 'rank_p2': r2['rank_p'], 'rank_aug_p2': r2['rank_aug_p'],
           'rank_Q_M': rM, 'rank_Q_Mb': rMb, 't_exact': tM + tMb,
           'rigorous_inconsistent_over_Q': bool(rMb == rM + 1)}
    print(json.dumps(rec), flush=True); log(rec)
