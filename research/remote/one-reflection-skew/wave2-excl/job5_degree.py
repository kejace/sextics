"""Y2-free (strict 5-variable) exclusion driver.  Usage: python3 job5_degree.py <d> <case1,case2,...> <msolve_threads> [tag]
cases: all | tor:XXXX (2-torsion class of the Y2-free sector) | sec:(a.b.c) | idx:i.j.k
Pipeline per case: exact (J_f:A)_d family with exact lifts -> mod-J_f Hessian quadrics + certified square killing ->
strict relaxed duals (5-var Koszul, no F) + certified square killing -> Macaulay / msolve zero-dim test (mod p => V_C={0}).
Writes d<d>_5v_fam.pkl, d<d>_5v_<case>.pkl/.ms/.lead and log lines (exp='5v-d<d>-<case>')."""
import sys, os, json, time, pickle, subprocess, re, math
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x, w2x5
from w2x import ctx
P = 2147483647; PMS = 1073741827
d = int(sys.argv[1]); cases = sys.argv[2].split(','); threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
tag = sys.argv[4] if len(sys.argv) > 4 else ''
t0 = time.time()
famf = f'd{d}_5v_fam.pkl'
if os.path.exists(famf):
    obj = pickle.load(open(famf, 'rb'))
    items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
    fam = w2x5.Family5(d, items)
else:
    fam = w2x5.Family5.from_e1(d, verbose=True)
    pickle.dump({'d': d, 'items': [(nm, D.to_dict(), [n.to_dict() for n in N], den, chi, kind) for (nm, D, N, den, chi, kind) in fam.items]}, open(famf, 'wb'))
    per = {}
    for it in fam.items: per[str(it[4])] = per.get(str(it[4]), 0) + 1
    w2x.logexp(exp=f'5v-d{d}-family', d=d, prime=None, family='(J_f:A)_d exact basis with exact lifts (Y2-free)', result=f'n={fam.n} J={sum(1 for it in fam.items if it[5]=="J")} extra={sum(1 for it in fam.items if it[5]=="extra")}',
               rank=fam.n, seconds=round(time.time()-t0, 1), job=f'job5_degree {d}', note=f'per sector {per}')
print('family', fam.n, round(time.time()-t0, 1), flush=True)
tab = {}
for chi in w2x.ALLCHR:
    ms = w2x5.secmons5(d, chi)
    if ms:
        e = ms[0]; tab[chi] = ''.join(str(x) for x in ((e[0]+e[1]+e[2])%2, (e[0]+e[2]+e[3])%2, (e[0]+e[2]+e[4])%2, e[5]%2))
def alive_for(case):
    if case == 'all': return set(range(fam.n))
    if case.startswith('tor:'): return {a for a in range(fam.n) if tab[fam.chi[a]] == case[4:]}
    if case.startswith('sec:'):
        chi = tuple(int(x) for x in case[4:].strip('()').split('.')); return {a for a in range(fam.n) if fam.chi[a] == chi}
    if case.startswith('idx:'): return {int(x) for x in case[4:].split('.')}
    raise ValueError(case)
def parse_leading(path):
    txt = open(path).read(); m = re.search(r'\[(.*)\]', txt, re.S)
    return [s.strip() for s in m.group(1).split(',') if s.strip()] if m else None
def zero_dim(mons, names):
    pure = {mo.split('^')[0] for mo in mons if '*' not in mo}
    return all(v in pure for v in names), sorted(set(names) - pure)

bankf = f'd{d}_5v_bank.pkl'
if os.path.exists(bankf): bank = pickle.load(open(bankf, 'rb'))
else:
    bank = w2x5.build_bank5(fam); pickle.dump(bank, open(bankf, 'wb'))
print('bank ready', round(time.time()-t0, 1), flush=True)
jd = w2x5.jac_duals_for5(fam)
for case in cases:
    t1 = time.time(); alive = alive_for(case)
    safe = case.replace(':', '_').replace('(', '').replace(')', '').replace(',', '.')
    rec = {'exp': f'5v-d{d}-{case}', 'd': d, 'case': case, 'alive0': len(alive), 'prime': P}
    if not alive:
        rec['result'] = 'empty'; w2x.logexp(**rec); print(rec, flush=True); continue
    prs = [(a, b) for a in sorted(alive) for b in sorted(alive) if a <= b]
    qJ = w2x5.quadrics_from_duals5(fam, jd, which_pairs=prs)
    alive1, hist1, an1 = w2x.kill_squares(qJ, alive, verbose=True, p=P)
    rec['modJ'] = {'hist': hist1, 'alive': len(alive1), 'nquads': len(qJ)}
    if not alive1:
        rec.update(result='EXCLUDED (mod-J_f squares)', seconds=round(time.time()-t1, 1)); w2x.logexp(**rec); print(json.dumps(rec), flush=True); continue
    rd = w2x5.relaxed_duals5(fam, subset=alive1, bank=bank, verbose=True)
    prs1 = [(a, b) for a in sorted(alive1) for b in sorted(alive1) if a <= b]
    qR = w2x5.quadrics_from_duals5(fam, rd, which_pairs=prs1)
    alive2, hist2, an2 = w2x.kill_squares(qR, alive1, verbose=True, p=P)
    rec['relaxed'] = {'hist': hist2, 'alive': len(alive2), 'nquads': len(qR), 'rank': an2['rank'], 'monos': len(an2['monos'])}
    pickle.dump({'alive': sorted(alive2), 'quads': qR, 'alive1': sorted(alive1), 'quadsJ': qJ}, open(f'd{d}_5v_{safe}{tag}.pkl', 'wb'))
    if not alive2:
        rec.update(result='EXCLUDED (relaxed squares)', seconds=round(time.time()-t1, 1)); w2x.logexp(**rec); print(json.dumps(rec), flush=True); continue
    s = len(alive2)
    if math.comb(s + 2, 3) <= 20000:
        r3, n3, _ = w2x.macaulay_full_rank(qR, alive2, 3, P); rec['macaulay3'] = [r3, n3]
        if r3 == n3:
            rec.update(result='EXCLUDED (Macaulay deg 3 full rank mod p)', seconds=round(time.time()-t1, 1)); w2x.logexp(**rec); print(json.dumps(rec), flush=True); continue
    msin = f'd{d}_5v_{safe}{tag}.ms'; msout = f'd{d}_5v_{safe}{tag}.lead'
    nl, nv = w2x.quads_to_msolve(qR, alive2, msin, p=PMS)
    rec['msolve'] = {'nquads': nl, 'nvars': nv}
    t2 = time.time()
    try:
        subprocess.run(['msolve', '-t', str(threads), '-g', '1', '-f', msin, '-o', msout], capture_output=True, text=True, timeout=8 * 3600)
        mons = parse_leading(msout); names = [f'l{a}' for a in sorted(alive2)]
        zd, missing = zero_dim(mons, names) if mons else (False, names)
        rec['msolve'].update(seconds=round(time.time()-t2, 1), n_leading=len(mons) if mons else None, zero_dim=zd, missing_pure_powers=len(missing))
        rec['result'] = 'EXCLUDED (strict relaxed ideal zero-dimensional mod p => V_C={0})' if zd else f'OPEN: positive-dimensional ({len(missing)} vars without pure power)'
        if not zd: rec['missing_vars'] = missing[:60]
    except subprocess.TimeoutExpired:
        rec['result'] = 'msolve timeout'
    rec['seconds'] = round(time.time()-t1, 1); rec['prime_msolve'] = PMS
    w2x.logexp(**rec); print(json.dumps(rec), flush=True)
print('DONE', round(time.time()-t0, 1), flush=True)
