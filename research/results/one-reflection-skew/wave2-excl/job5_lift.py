"""Lifted (second-level) strict relaxation on a Y2-free sub-family.  Usage: python3 job5_lift.py <d> <case> <threads>
Loads d<d>_5v_fam.pkl and d<d>_5v_<case>.pkl (alive2 + first-level quadrics), builds the K-indexed bank, computes the
lifted cubics on the surviving alive set, writes msolve input (quadrics + cubics) and tests zero-dimensionality."""
import sys, os, json, time, pickle, subprocess, re
sys.path.insert(0, '/tmp/hodge/results/one-reflection-skew/wave2-excl')
import w2x, w2x5
from w2x import ctx
PMS = 1073741827
d = int(sys.argv[1]); case = sys.argv[2]; threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
safe = case.replace(':', '_').replace('(', '').replace(')', '').replace(',', '.')
obj = pickle.load(open(f'd{d}_5v_fam.pkl', 'rb'))
items = [(nm, ctx.from_dict(D), [ctx.from_dict(n) for n in N], den, tuple(chi), kind) for (nm, D, N, den, chi, kind) in obj['items']]
fam = w2x5.Family5(d, items)
leaf = pickle.load(open(f'd{d}_5v_{safe}.pkl', 'rb'))
alive = set(leaf['alive']); quads = leaf['quads']
t0 = time.time()
bkf = f'd{d}_5v_bankK.pkl'
if os.path.exists(bkf): bank, Ks = pickle.load(open(bkf, 'rb'))
else:
    bank, Ks = w2x5.build_bank5k(fam); pickle.dump((bank, Ks), open(bkf, 'wb'))
print('bankK ready', round(time.time()-t0, 1), flush=True)
cubics = w2x5.lift2_cubics5(fam, alive, bank, Ks, verbose=True)
print('cubics', len(cubics), round(time.time()-t0, 1), flush=True)
pickle.dump({'alive': sorted(alive), 'cubics': cubics}, open(f'd{d}_5v_{safe}_lift.pkl', 'wb'))
msin = f'd{d}_5v_{safe}_lift.ms'; msout = f'd{d}_5v_{safe}_lift.lead'
nl, nv = w2x5.polys_to_msolve(quads, cubics, alive, msin, p=PMS)
t1 = time.time()
subprocess.run(['msolve', '-t', str(threads), '-g', '1', '-f', msin, '-o', msout], capture_output=True, text=True, timeout=6 * 3600)
txt = open(msout).read(); m = re.search(r'\[(.*)\]', txt, re.S)
mons = [s.strip() for s in m.group(1).split(',') if s.strip()] if m else []
pure = {mo.split('^')[0] for mo in mons if '*' not in mo}
names = [f'l{a}' for a in sorted(alive)]
missing = sorted(set(names) - pure)
rec = {'exp': f'5v-d{d}-{case}-lift2', 'd': d, 'case': case, 'alive': len(alive), 'nquads': len(w2x.restrict_quads(quads, alive)), 'ncubics': len(cubics),
       'npolys_msolve': nl, 'n_leading': len(mons), 'zero_dim': not missing, 'missing_vars': missing[:40], 'prime_msolve': PMS,
       'result': 'EXCLUDED (lifted strict relaxed ideal zero-dimensional mod p => V_C={0})' if not missing else f'OPEN after lift2 ({len(missing)} vars without pure power)',
       'seconds': round(time.time()-t0, 1), 'msolve_seconds': round(time.time()-t1, 1)}
w2x.logexp(**rec); print(json.dumps(rec), flush=True)
