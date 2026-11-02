"""Served-target spaces (target_image mod P2) for structured irreducible-over-Q blocks of degree <= dmax; span test for A."""
import sys, time, json
from w2lib import *
dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
f0 = add(power(S0, 3), scale(mul(mul(P0, P0), S0), 2))
rad = [add(mul(X(i), FX[i]), mul(X(i + 3), FX[i + 3])) for i in range(3)]
T0 = add(power(X0, 4), power(Y0, 4))
blocks = [('X0', X0), ('Y0', Y0), ('X1', X1), ('Y1', Y1), ('X2', X2), ('Y2', Y2),
          ('S0', S0), ('delta', DELTA), ('S1', S1), ('S2', S2), ('S0+S2', add(S0, S2)), ('S0-S2', sub(S0, S2)), ('X0^2+X2^2', add(mul(X0, X0), mul(X2, X2))),
          ('S0+S1', add(S0, S1)), ('X0^2-X1^2', sub(mul(X0, X0), mul(X1, X1))), ('X1^2-Y1^2', sub(mul(X1, X1), mul(Y1, Y1))), ('X2^2-Y2^2', sub(mul(X2, X2), mul(Y2, Y2))),
          ('GS', GS), ('GP', GP), ('Q1', Q1), ('T0', T0), ('S0^2+2P0^2', add(mul(S0, S0), scale(mul(P0, P0), 2))), ('Q2', Q2), ('X1^4+Y1^4', add(power(X1, 4), power(Y1, 4))),
          ('F_X0', FX[0]), ('F_Y0', FX[3]), ('F_X1', FX[1]), ('F_Y1', FX[4]), ('F_X2', FX[2]), ('F_Y2', FX[5]),
          ('f0', f0), ('rad0', rad[0]), ('rad1', rad[1]), ('rad2', rad[2]), ('X1^6+Y1^6', add(power(X1, 6), power(Y1, 6))), ('X2^6+Y2^6', add(power(X2, 6), power(Y2, 6))),
          ('D8*', D8)]
served_file = PK + '/served_local.json'
try: served = json.load(open(served_file))
except Exception: served = {}
for nm, D in blocks:
    if core.degree(D) > dmax or nm in served: continue
    t0 = time.time()
    if core.parity_of(D) is None:
        log('served', D=nm, d=core.degree(D), result='skipped-mixed-parity'); continue
    rank, basis, S = fixed_d.target_image(F, q, D, p=p, verbose=False)
    served[nm] = {'d': core.degree(D), 'basis': [core.to_str(B) for B in basis]}
    json.dump(served, open(served_file, 'w'))
    amons = sorted({e for B in basis for e in B} | set(Ap))
    inA = False
    if basis:
        M = nmod_mat(len(basis), len(amons), [B.get(e, 0) % p for B in basis for e in amons], p)
        M2 = nmod_mat(len(basis) + 1, len(amons), [B.get(e, 0) % p for B in basis + [Ap] for e in amons], p)
        inA = (M.rank() == M2.rank())
    log('served', D=nm, d=core.degree(D), dim_served=rank, A_in_served=inA, rows=S.nrows, cols=S.ncols, seconds=round(time.time() - t0, 1),
        note='target_image mod P2 (local): space of Y-odd degree-6 targets served by this denominator at q=1')
# span test over everything computed so far
rows = []
for nm, v in served.items():
    rows.extend(core.parse(s, mod=p) if isinstance(s, str) else s for s in v['basis'])
amons = sorted({e for B in rows for e in B} | set(Ap))
if rows:
    M = nmod_mat(len(rows), len(amons), [B.get(e, 0) % p for B in rows for e in amons], p)
    M2 = nmod_mat(len(rows) + 1, len(amons), [B.get(e, 0) % p for B in rows + [Ap] for e in amons], p)
    log('served-span', blocks=sorted(served.keys()), dim_span=M.rank(), A_in_span=(M.rank() == M2.rank()),
        note='sum principle: served(prod D_i) = sum served(D_i) (support depends only on the irreducible factors); A in span <=> some product of the listed blocks could serve A')
