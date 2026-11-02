"""Degree-d exclusion driver (exact):  python3 job_degree.py <d> [stages]
A: E1 family (exact kernel basis, exact lifts)  B: mod-J Hessian quadrics + certified square killing
C: relaxed mod-F quadrics (all-character Koszul corrections) on survivors + certified square killing
D: leaf analysis: Macaulay rank (small) / msolve+Singular inputs.  Checkpoints in wave2-excl/d<d>_*."""
import sys, os, json, time, pickle
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx

d = int(sys.argv[1]); stages = sys.argv[2] if len(sys.argv) > 2 else 'ABCD'
P = 2147483647
pre = f'd{d}'
def save(name, obj):
    with open(f'{pre}_{name}.pkl', 'wb') as fh: pickle.dump(obj, fh)
def load(name):
    with open(f'{pre}_{name}.pkl', 'rb') as fh: return pickle.load(fh)
def fam_to_pickle(fam):
    return {'d': fam.d, 'items': [(nm, D.to_dict(), [n.to_dict() for n in N], den, chi, kind) for (nm, D, N, den, chi, kind) in fam.items]}
def fam_from_pickle(obj):
    items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
    return w2x.Family(obj['d'], items)

t0 = time.time()
if 'A' in stages:
    fam = w2x.Family.from_e1(d, gauge=True, verbose=True)
    save('fam', fam_to_pickle(fam))
    kinds = {k: sum(1 for it in fam.items if it[5] == k) for k in ('J', 'extra')}
    persec = {}
    for it in fam.items: persec[str(it[4])] = persec.get(str(it[4]), 0) + 1
    w2x.logexp(exp=f'stageA-family-d{d}', family='gauge-fixed exact E1 kernel basis with exact lifts', d=d, prime=None,
               result=f'n={fam.n} J={kinds["J"]} extra={kinds["extra"]}', rank=fam.n, seconds=round(time.time()-t0, 1), job=f'job_degree {d}',
               note=f'per-sector counts {persec}; lifts verified exactly (sum N_i F_i = den A D)')
    print('A done', fam.n, kinds, round(time.time()-t0, 1), flush=True)
else:
    fam = fam_from_pickle(load('fam'))
if 'B' in stages:
    t1 = time.time()
    duals = w2x.jac_duals_for(fam)
    quads = w2x.quadrics_from_duals(fam, duals)
    save('quads_modJ', quads)
    alive, hist, an = w2x.kill_squares(quads, range(fam.n), verbose=True, p=P)
    res = {'alive': sorted(alive), 'hist': hist, 'monos': an['monos'], 'rank_final': an['rank'],
           'alive_names': [fam.items[a][0] for a in sorted(alive)], 'alive_kinds': [fam.items[a][5] for a in sorted(alive)],
           'alive_sectors': sorted({str(fam.items[a][4]) for a in alive})}
    with open(f'{pre}_stageB.json', 'w') as fh: json.dump(res, fh)
    w2x.logexp(exp=f'stageB-modJ-d{d}', family='Hessian relaxation mod J on the gauge-fixed E1 kernel', d=d, prime=P,
               result=f'{fam.n} -> {len(alive)} survivors after certified square killing; {len(an["monos"])} off-diagonal monomials in span', rank=an['rank'],
               seconds=round(time.time()-t1, 1), job=f'job_degree {d}', note=f'history {hist}; squares certified exactly over Q; sectors of survivors {res["alive_sectors"]}')
    print('B done', len(alive), hist, round(time.time()-t0, 1), flush=True)
else:
    quads = load('quads_modJ'); res = json.load(open(f'{pre}_stageB.json')); alive = set(res['alive'])
if 'C' in stages and alive:
    t1 = time.time()
    rd = w2x.relaxed_duals(fam, subset=alive, verbose=True)
    prs = [(a, b) for a in sorted(alive) for b in sorted(alive) if a <= b]
    quads2 = w2x.quadrics_from_duals(fam, rd, which_pairs=prs)
    save('quads_relaxed', quads2)
    alive2, hist2, an2 = w2x.kill_squares(quads2, alive, verbose=True, p=P)
    res2 = {'alive': sorted(alive2), 'hist': hist2, 'monos': an2['monos'], 'rank_final': an2['rank'],
            'alive_names': [fam.items[a][0] for a in sorted(alive2)], 'alive_kinds': [fam.items[a][5] for a in sorted(alive2)],
            'alive_sectors': sorted({str(fam.items[a][4]) for a in alive2})}
    with open(f'{pre}_stageC.json', 'w') as fh: json.dump(res2, fh)
    w2x.logexp(exp=f'stageC-relaxed-d{d}', family='relaxed mod-F system (exact duals, all-character Koszul corrections) on mod-J survivors', d=d, prime=P,
               result=f'{len(alive)} -> {len(alive2)} survivors; {len(an2["monos"])} off-diagonal monomials in span; quadric rank {an2["rank"]}', rank=an2['rank'],
               seconds=round(time.time()-t1, 1), job=f'job_degree {d}', note=f'history {hist2}; survivors sectors {res2["alive_sectors"]}')
    print('C done', len(alive2), hist2, round(time.time()-t0, 1), flush=True)
    alive = alive2
elif 'C' in stages:
    quads2 = []
else:
    quads2 = load('quads_relaxed'); res2 = json.load(open(f'{pre}_stageC.json')); alive = set(res2['alive'])
if 'D' in stages and alive:
    t1 = time.time()
    s = len(alive)
    import math
    nc3 = math.comb(s + 2, 3)
    if nc3 <= 30000:
        r3, n3, rows3 = w2x.macaulay_full_rank(quads2, alive, 3, P)
        note = f'Macaulay deg3 rank {r3}/{n3} (rows {rows3})'
        if r3 < n3 and math.comb(s + 3, 4) <= 60000:
            r4, n4, rows4 = w2x.macaulay_full_rank(quads2, alive, 4, P)
            note += f'; deg4 {r4}/{n4}'
            full = (r4 == n4)
        else:
            full = (r3 == n3)
        w2x.logexp(exp=f'stageD-macaulay-d{d}', family='leaf: relaxed quadrics on survivors', d=d, prime=P, result='V={0} (full column rank) => EXCLUDED' if full else 'not full rank; needs Groebner/case split',
                   rank=r3, seconds=round(time.time()-t1, 1), job=f'job_degree {d}', note=note)
        print('D', note, 'EXCLUDED' if full else 'OPEN', flush=True)
    nl, nv = w2x.quads_to_msolve(quads2, alive, f'{pre}_leaf.ms', p=1073741827)
    # Singular input too
    al = sorted(alive)
    with open(f'{pre}_leaf.sing', 'w') as fh:
        fh.write('ring r = 32003, (' + ','.join(f'l{a}' for a in al) + '), dp;\n')
        polys = []
        for qd in w2x.restrict_quads(quads2, al):
            terms = []
            for (a, b), v in qd['coef'].items():
                vv = v % 32003
                if vv: terms.append(f'{vv}*l{a}*l{b}')
            if terms: polys.append('+'.join(terms))
        fh.write('ideal I = ' + ',\n'.join(polys) + ';\n')
        fh.write('option(redSB); ideal G = std(I); print("DIM"); print(dim(G)); print("HILB"); hilb(G); print("SIZE"); print(size(G)); quit;\n')
    print('D wrote leaf inputs', nl, 'quadrics', nv, 'vars', flush=True)
print('ALL DONE', round(time.time()-t0, 1), flush=True)
