"""Job: excl-fixedD.  (1) sector fixed-D re-solves at P2 for D8*, Z_i*D8* (pure characters); (2) engine cross-check
(fixed_d.solve_fixed_D_fast, different method) at P2 for D8*, X1*D8*; (3) served-target spaces (target_image) for irreducible
structured blocks and the span test A in sum served(D_i) (sum principle)."""
import sys, time, json
from w2lib import *
OUT = PK + '/fixedD'
os.makedirs(OUT, exist_ok=True)
parts = sys.argv[1:] or ['1', '2', '3']

if '1' in parts:
    test_D(D8, 'D8*', 'wave1-replay', note='fqlx/wave1: inconsistent (rank 718, nullity 109 at p=2^31-1)')
    for nm, Z in (('X1', X1), ('Y1', Y1), ('X2', X2), ('Y2', Y2), ('X0', X0), ('Y0', Y0)):
        test_D(mul(D8, Z), f'{nm}*D8*', 'wave1-replay', note='wave1: inconsistent at p=2^31-1')

if '2' in parts:
    for nm, D in (('D8*', D8), ('X1*D8*', mul(D8, X1))):
        t0 = time.time()
        ok, info = fixed_d.e1_solvable(F, A, q, D, p=p)
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, p=p, verbose=True, want_exact=False)
        log('fixedD-engine', D=nm, family='engine-crosscheck', d=core.degree(D), e1=ok, result='CONSISTENT' if out['consistent'] else 'inconsistent',
            rank=out['rank'], nullity=out['nullity'], rows=out['rows'], cols=out['cols'], seconds=round(time.time() - t0, 1),
            note='hodge.fixed_d.solve_fixed_D_fast (compressed modular rref, parity-split, no character restriction)')

if '3' in parts:
    f0 = add(power(S0, 3), scale(mul(mul(P0, P0), S0), 2))
    rad = [add(mul(X(i), FX[i]), mul(X(i + 3), FX[i + 3])) for i in range(3)]
    T0 = add(power(X0, 4), power(Y0, 4))
    blocks = [('X0', X0), ('Y0', Y0), ('X1', X1), ('Y1', Y1), ('X2', X2), ('Y2', Y2),
              ('S0', S0), ('delta', DELTA), ('S1', S1), ('S2', S2), ('S0+S2', add(S0, S2)), ('S0-S2', sub(S0, S2)), ('X0^2+X2^2', add(mul(X0, X0), mul(X2, X2))),
              ('S0+S1', add(S0, S1)), ('X0^2-X1^2', sub(mul(X0, X0), mul(X1, X1))),
              ('GS', GS), ('GP', GP), ('Q1', Q1), ('T0', T0), ('S0^2+2P0^2', add(mul(S0, S0), scale(mul(P0, P0), 2))), ('Q2', Q2),
              ('f0', f0), ('rad0', rad[0]), ('rad1', rad[1]), ('rad2', rad[2]), ('X1^6+Y1^6', add(power(X1, 6), power(Y1, 6))), ('X2^6+Y2^6', add(power(X2, 6), power(Y2, 6))),
              ('F_X0', FX[0]), ('F_Y0', FX[3]), ('F_X1', FX[1]), ('F_X2', FX[2]),
              ('D8*', D8)]
    served = {}
    span_rows = []
    for nm, D in blocks:
        t0 = time.time()
        if core.parity_of(D) is None:
            log('served', D=nm, d=core.degree(D), result='skipped-mixed-parity'); continue
        try:
            rank, basis, S = fixed_d.target_image(F, q, D, p=p, verbose=False)
        except Exception as ex:
            log('served', D=nm, d=core.degree(D), result='ERROR ' + repr(ex)); continue
        served[nm] = [core.to_str(B) for B in basis]
        json.dump(served, open(OUT + '/served.json', 'w'))
        # is A in served(D)?  and accumulate the span
        amons = sorted({e for B in basis for e in B} | set(Ap))
        inA = False
        if basis:
            aidx = {e: i for i, e in enumerate(amons)}
            M = nmod_mat(len(basis), len(amons), [B.get(e, 0) % p for B in basis for e in amons], p)
            M2 = nmod_mat(len(basis) + 1, len(amons), [B.get(e, 0) % p for B in basis + [Ap] for e in amons], p)
            inA = (M.rank() == M2.rank())
        span_rows.extend(basis)
        log('served', D=nm, d=core.degree(D), dim_served=rank, A_in_served=inA, rows=S.nrows, cols=S.ncols, seconds=round(time.time() - t0, 1),
            note='target_image mod p: space of Y-odd degree-6 targets served by this denominator (q=1)')
    # total span
    amons = sorted({e for B in span_rows for e in B} | set(Ap))
    if span_rows:
        M = nmod_mat(len(span_rows), len(amons), [B.get(e, 0) % p for B in span_rows for e in amons], p)
        M2 = nmod_mat(len(span_rows) + 1, len(amons), [B.get(e, 0) % p for B in span_rows + [Ap] for e in amons], p)
        log('served-span', blocks=[nm for nm, _ in blocks if nm in served], dim_span=M.rank(), A_in_span=(M.rank() == M2.rank()),
            note='sum principle: served(prod D_i) = sum served(D_i); A in span <=> some product of the listed blocks could serve A')
print('DONE', file=sys.stderr)
