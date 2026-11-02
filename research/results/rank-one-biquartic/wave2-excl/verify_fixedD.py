"""V3: re-verify fixed-D inconsistency of the structural denominators with the DENSE solver at a NEW prime
(FixedDSystem.solve = full rref, no row compression).  Parallel over denominators."""
import sys, json, time, multiprocessing as mp
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
p = int(sys.argv[1]) if len(sys.argv) > 1 else 2147483423
LOGF = f'{PK}/verify_fixedD_{p}.log.jsonl'
DX = {}; DY = {}
for i in range(3): DX = add(DX, mul(EULER[i], Fi[i]))
for i in range(3, 6): DY = add(DY, mul(EULER[i], Fi[i]))
cands = [('E_X(F)', DX), ('E_Y(F)', DY), ('E_X(F)-E_Y(F)', sub(DX, DY)), ('X0*F_0', mul(EULER[0], Fi[0])), ('Y0*F_3', mul(EULER[3], Fi[3])),
         ('X1*F_1', mul(EULER[1], Fi[1])), ('X0*F_1', mul(EULER[0], Fi[1])), ('Y0*F_4', mul(EULER[3], Fi[4])),
         ('F_X0', Fi[0]), ('F_X1', Fi[1]), ('F_X2', Fi[2]), ('F_X0+F_X1+F_X2', add(add(Fi[0], Fi[1]), Fi[2])),
         ('F_Y0', Fi[3]), ('F_Y1', Fi[4]), ('F_Y2', Fi[5]),
         ('f', f), ('g', g), ('f+g', add(f, g)), ('f-g', sub(f, g)), ('L', L), ('M', M), ('L^2', mul(L, L)), ('L*M', mul(L, M)), ('M^2', mul(M, M)),
         ('L^2*M', L2M), ('L*M^2', mul(mul(L, M), M)), ('L^3', mul(mul(L, L), L)), ('M^3', mul(mul(M, M), M))]
def run(item):
    name, D = item
    t0 = time.time()
    try:
        ok, info = fixed_d.e1_solvable(F, A, q, D, p)
        rec_ = dict(exp='reverify-fixedD-dense', D=name, d=core.degree(D), prime=p, e1=bool(ok), e1_rank=info['rank'], e1_rank_aug=info['rank_aug'])
        if ok:
            S = fixed_d.FixedDSystem(F, A, q, D, p)
            r = S.solve()
            rec_.update(result='consistent' if r['consistent'] else 'inconsistent', rank=r['rank'], nullity=r['nullity'], rows=S.nrows, cols=S.ncols)
        else:
            rec_.update(result='E1 fails (A*D not in J mod p)')
    except Exception as ex:
        rec_ = dict(exp='reverify-fixedD-dense', D=name, prime=p, result=f'error {ex!r}')
    rec_['seconds'] = round(time.time() - t0, 1)
    rec_['note'] = 'dense rref of the full cleared system (no compression) at a new prime'
    return rec_
if __name__ == '__main__':
    with mp.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 8) as pool:
        for rec_ in pool.imap_unordered(run, cands):
            print(json.dumps(rec_), flush=True); log(rec_, path=LOGF)
    print('done')
