"""Degree-7 sub-family exclusion via relaxed quadrics + msolve Groebner (zero-dimensionality => V_C = {0}).
Usage: python3 job_d7_sub.py <d> <case1,case2,...> <threads>
cases: odd (all Y2-odd sectors, mixtures), even (all Y2-even sectors, mixtures), tor:XXXX (pure 2-torsion class),
       full (all sectors), sec:(a,b,c) (single fine sector).  Each case: own relaxed duals on its alive set."""
import sys, os, json, time, pickle, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x
from w2x import ctx
P = 2147483647; PMS = 1073741827
d = int(sys.argv[1]); cases = sys.argv[2].split(','); threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
tag = sys.argv[4] if len(sys.argv) > 4 else ''
obj = pickle.load(open(f'd{d}_fam.pkl', 'rb'))
items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
fam = w2x.Family(d, items)
tab = {}
for chi in w2x.ALLCHR:
    ms = w2x.secmons(d, chi)
    if not ms: continue
    e = ms[0]
    tab[chi] = {'e5': e[5] % 2, 'tor': ''.join(str(x) for x in ((e[0]+e[1]+e[2])%2, (e[0]+e[2]+e[3])%2, (e[0]+e[2]+e[4])%2, e[5]%2))}

def alive_for(case):
    if case == 'full': return set(range(fam.n))
    if case == 'odd': return {a for a in range(fam.n) if tab[fam.chi[a]]['e5'] == 1}
    if case == 'even': return {a for a in range(fam.n) if tab[fam.chi[a]]['e5'] == 0}
    if case.startswith('tor:'): return {a for a in range(fam.n) if tab[fam.chi[a]]['tor'] == case[4:]}
    if case.startswith('sec:'):
        chi = tuple(int(x) for x in case[4:].strip('()').split('.'))
        return {a for a in range(fam.n) if fam.chi[a] == chi}
    if case.startswith('idx:'): return {int(x) for x in case[4:].split('.')}
    raise ValueError(case)

def parse_leading(path):
    txt = open(path).read()
    m = re.search(r'\[(.*)\]', txt, re.S)
    if not m: return None
    mons = [s.strip() for s in m.group(1).split(',') if s.strip()]
    return mons

def zero_dim(mons, names):
    pure = set()
    for mo in mons:
        parts = mo.split('*')
        if len(parts) == 1:
            v = parts[0].split('^')[0]
            pure.add(v)
    return all(v in pure for v in names), sorted(set(names) - pure)

for case in cases:
    t0 = time.time()
    alive = alive_for(case)
    rec = {'exp': f'd{d}-subfamily-{case}', 'd': d, 'case': case, 'alive0': len(alive)}
    if not alive:
        rec.update(result='empty family'); w2x.logexp(**rec); print(rec, flush=True); continue
    # mod-J quadrics restricted (reuse full ones if present)
    if os.path.exists(f'd{d}_quads_modJ.pkl'):
        qJ = w2x.restrict_quads(pickle.load(open(f'd{d}_quads_modJ.pkl', 'rb')), alive)
    else:
        qJ = w2x.quadrics_from_duals(fam, w2x.jac_duals_for(fam), which_pairs=[(a, b) for a in sorted(alive) for b in sorted(alive) if a <= b])
    alive1, hist1, an1 = w2x.kill_squares(qJ, alive, verbose=True, p=P)
    rec['modJ'] = {'hist': hist1, 'alive': len(alive1)}
    if not alive1:
        rec.update(result='EXCLUDED by mod-J squares', seconds=round(time.time()-t0, 1)); w2x.logexp(**rec); print(rec, flush=True); continue
    rd = w2x.relaxed_duals(fam, subset=alive1, verbose=False)
    prs = [(a, b) for a in sorted(alive1) for b in sorted(alive1) if a <= b]
    qR = w2x.quadrics_from_duals(fam, rd, which_pairs=prs)
    alive2, hist2, an2 = w2x.kill_squares(qR, alive1, verbose=True, p=P)
    rec['relaxed'] = {'hist': hist2, 'alive': len(alive2), 'nquads': len(qR), 'rank': an2['rank'], 'monos': len(an2['monos'])}
    pickle.dump({'alive': sorted(alive2), 'quads': qR}, open(f'd{d}_sub_{case.replace(":", "_").replace("(", "").replace(")", "").replace(",", ".")}{tag}.pkl', 'wb'))
    if not alive2:
        rec.update(result='EXCLUDED by relaxed squares', seconds=round(time.time()-t0, 1)); w2x.logexp(**rec); print(rec, flush=True); continue
    s = len(alive2)
    import math
    if math.comb(s + 2, 3) <= 20000:
        r3, n3, _ = w2x.macaulay_full_rank(qR, alive2, 3, P)
        rec['macaulay3'] = [r3, n3]
        if r3 == n3:
            rec.update(result='EXCLUDED (Macaulay deg 3 full rank)', seconds=round(time.time()-t0, 1)); w2x.logexp(**rec); print(rec, flush=True); continue
    safe = case.replace(':', '_').replace('(', '').replace(')', '').replace(',', '.')
    msin = f'd{d}_sub_{safe}{tag}.ms'; msout = f'd{d}_sub_{safe}{tag}.lead'
    nl, nv = w2x.quads_to_msolve(qR, alive2, msin, p=PMS)
    rec['msolve'] = {'nquads': nl, 'nvars': nv}
    t1 = time.time()
    try:
        cp = subprocess.run(['msolve', '-t', str(threads), '-g', '1', '-f', msin, '-o', msout], capture_output=True, text=True, timeout=6 * 3600)
        mons = parse_leading(msout)
        names = [f'l{a}' for a in sorted(alive2)]
        zd, missing = zero_dim(mons, names) if mons else (False, names)
        rec['msolve'].update(seconds=round(time.time()-t1, 1), n_leading=len(mons) if mons else None, zero_dim=zd, missing_pure_powers=len(missing))
        rec['result'] = 'EXCLUDED (relaxed ideal zero-dimensional mod p => V_C={0})' if zd else f'OPEN: relaxed ideal positive-dimensional ({len(missing)} vars without pure power)'
        if not zd: rec['missing_vars'] = missing[:40]
    except subprocess.TimeoutExpired:
        rec['result'] = 'msolve timeout'
    rec['seconds'] = round(time.time()-t0, 1)
    rec['prime'] = PMS
    w2x.logexp(**rec); print(json.dumps(rec), flush=True)
print('DONE', flush=True)
