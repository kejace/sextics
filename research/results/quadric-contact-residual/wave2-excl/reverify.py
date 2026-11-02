"""Job: excl-reverify.  Adversarial re-verification of wave-1 / fqlx claims at a new prime (P2), own code.
Parts: (1) E1 colon dims d=0..9 per sector; (2) d=5 Hessian columns rank; (3) d=6 full closure + Singular;
(4) d=7 parity sectors (8) closure + Singular; (5) all pure G_aff sectors at d=8 and d=9 closure + Singular."""
import sys, time, json
from w2lib import *
from family import *
OUT = PK + '/reverify'
os.makedirs(OUT, exist_ok=True)
parts = sys.argv[1:] or ['1', '2', '3', '4', '5']

if '1' in parts:
    known = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 6, 6: 36, 7: 126, 8: 337, 9: 760}
    for d in range(0, 10):
        t0 = time.time()
        prof, extras = e1_profile(d, verbose=False)
        colon = sum(v['kernel'] for v in prof.values()); dimJ = sum(v['dimJ'] for v in prof.values())
        extra = {k: [core.to_str(P)[:80] for P in v] for k, v in extras.items()}
        log('E1-colon', d=d, dim_colon=colon, dim_J=dimJ, n_extra=colon - dimJ, agrees_with_fqlx=(colon == known[d]),
            extra_sectors=list(extras.keys()), seconds=round(time.time() - t0, 1), note='dim (J:A)_d over all G_aff sectors; fqlx/dossier/wave1 values 0,0,0,0,0,6,36,126,337,760')
        json.dump({'profile': prof, 'extras': extra}, open(OUT + f'/e1_profile_d{d}.json', 'w'))

if '2' in parts:
    t0 = time.time()
    # d=5: D = sum c_i F_i, Hessian columns A*d_ij F in R_10[chi]
    Ds = [Fi[i] for i in range(NV)]
    Ns = []
    for i in range(NV):
        N = [{} for _ in range(NV)]; N[i] = dict(Ap); Ns.append(N)
    quads, qi = quadrics_from_family(Ds, Ns, 5, verbose=False)
    rank, forced, pairs = quad_rank_and_squares(quads, 6)
    neq, nv = write_systems(quads, 6, OUT + '/d5')
    v = singular_verdict(OUT + '/d5.sing')
    log('closure-d5', d=5, n_unknowns=6, n_quadrics=len(quads), quadric_rank=rank, n_pairs=len(pairs), forced_zero_by_squares=forced,
        singular_dim=v['sdim'], result='EXCLUDED' if v['sdim'] == 0 else 'OPEN', seconds=round(time.time() - t0, 1),
        note='fqlx: 21 Hessian columns have rank 21 in R_10 -> all c_i c_j = 0')

if '3' in parts:
    t0 = time.time()
    # d=6: D = sum m_ij Z_i F_j (36), gauge m_55 (D = Y2 F_Y2) dropped
    Ds, Ns, names = [], [], []
    for i in range(NV):
        for j in range(NV):
            Ds.append(mul({tuple(1 if t == i else 0 for t in range(NV)): 1}, Fi[j], p))
            N = [{} for _ in range(NV)]; N[j] = mul(Ap, {tuple(1 if t == i else 0 for t in range(NV)): 1}, p); Ns.append(N)
            names.append((i, j))
    quads, qi = quadrics_from_family(Ds, Ns, 6, verbose=False)
    gauge = names.index((5, 5))
    keep = [k for k in range(36) if k != gauge]
    kidx = {k: i for i, k in enumerate(keep)}
    quads2 = []
    for cf in quads:
        cf2 = {(kidx[k], kidx[l]): v for (k, l), v in cf.items() if k != gauge and l != gauge}
        if cf2: quads2.append(cf2)
    rank, forced, pairs = quad_rank_and_squares(quads2, 35)
    neq, nv = write_systems(quads2, 35, OUT + '/d6')
    v = singular_verdict(OUT + '/d6.sing')
    log('closure-d6', d=6, n_unknowns=35, n_quadrics=len(quads2), quadric_rank=rank, n_pairs=len(pairs), n_forced_by_squares=len(forced),
        forced_names=[names[keep[k]] for k in forced], singular_dim=v['sdim'], singular_vdim=v['svdim'],
        result='EXCLUDED' if v['sdim'] == 0 else 'OPEN', seconds=round(time.time() - t0, 1),
        note='fqlx: 630 columns rank 334 in R_12; 28 off-diagonal m_ij killed by squares; remaining 7 killed by explicit relations')

if '4' in parts:
    t0 = time.time()
    # d=7: D = sum c_{m,j} m F_j, m in S_2 (126); gauge (Z_k Z_5, j=5) dropped (6) -> 120; split by parity sector (b0,b2,y) of chi(D)
    Ds, Ns, names = [], [], []
    for m in monomials(2):
        for j in range(NV):
            if j == 5 and m[5] >= 1: continue
            Ds.append(shift(Fi[j], m, 1, p))
            N = [{} for _ in range(NV)]; N[j] = mul(Ap, {m: 1}, p); Ns.append(N); names.append((m, j))
    assert len(Ds) == 120
    cD = [chi_of(D) for D in Ds]
    par = [(c[0], c[1], c[3]) for c in cD]
    sectors = sorted(set(par))
    quads, qi = quadrics_from_family(Ds, Ns, 7, verbose=False)
    rank, forced, pairs = quad_rank_and_squares(quads, 120)
    neq, nv = write_systems(quads, 120, OUT + '/d7mixed')
    json.dump({'names': [[list(m), j] for m, j in names], 'chi': [list(c) for c in cD], 'quads': [[[k, l, v] for (k, l), v in cf.items()] for cf in quads]}, open(OUT + '/d7mixed.quads.json', 'w'))
    log('closure-d7-mixed-stats', d=7, n_unknowns=120, n_quadrics=len(quads), quadric_rank=rank, n_pairs=len(pairs), n_forced_by_squares=len(forced),
        seconds=round(time.time() - t0, 1), note='full 120-unknown system (all sectors mixed) exported to reverify/d7mixed.{ms,sing}; dossier: rank 950, no squares killed')
    for s in sectors:
        t1 = time.time()
        keep = [k for k in range(120) if par[k] == s]
        kidx = {k: i for i, k in enumerate(keep)}
        quads2 = []
        for cf in quads:
            cf2 = {(kidx[k], kidx[l]): v for (k, l), v in cf.items() if k in kidx and l in kidx}
            if cf2: quads2.append(cf2)
        r2, f2, p2 = quad_rank_and_squares(quads2, len(keep))
        tag = OUT + '/d7par_' + ''.join(map(str, s))
        neq, nv = write_systems(quads2, len(keep), tag)
        v = singular_verdict(tag + '.sing')
        log('closure-d7-parity-sector', d=7, parity_sector=str(s), n_unknowns=len(keep), n_quadrics=len(quads2), quadric_rank=r2, n_forced_by_squares=len(f2),
            singular_dim=v['sdim'], singular_vdim=v['svdim'], result='EXCLUDED' if v['sdim'] == 0 else 'OPEN', seconds=round(time.time() - t1, 1),
            note='fqlx: (1,1,*) sectors 12 params killed by squares; others 14 of 16 killed then u^2+v^2=uv=0')

if '5' in parts:
    for d in (8, 9):
        chars = sorted({chi_of({e: 1}) for e in monomials(d)})
        for c in chars:
            t0 = time.time()
            tag = OUT + f'/d{d}_chi' + ''.join(map(str, c))
            info = sector_closure(d, c, tag, run_singular=True, timeout=7200, verbose=False)
            log('closure-pure-sector', d=d, chi=str(c), **{k: v for k, v in info.items() if k not in ('d', 'chi')},
                note='gauge-fixed Hessian relaxation of the full E1-admissible family in this G_aff character sector; Singular dim 0 = only c=0')
print('DONE', file=sys.stderr)
