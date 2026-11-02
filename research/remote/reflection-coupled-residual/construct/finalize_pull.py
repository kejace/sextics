"""After the remote chains finish: pull the lane dir and append one log line per finished experiment to the lane log.
Usage: python3 finalize_pull.py   (run locally; uses rjob pull)"""
import json, os, subprocess, glob, time
LANE = 'reflection-coupled-residual'
subprocess.run(['python3', '-m', 'hodge.rjob', 'pull', '--lane', LANE], cwd='/home/kejace/devel/sextics-work')
R = f'/tmp/hodge/remote/{LANE}/construct'
L = f'/tmp/hodge/results/{LANE}/wave2-construct.log.jsonl'
done = set()
for line in open(L):
    try:
        d = json.loads(line); done.add(d.get('exp') + '|' + str(d.get('family')))
    except Exception: pass
out = []
for js in glob.glob(f'{R}/*.json'):
    name = os.path.basename(js)[:-5]
    try: S = json.load(open(js))
    except Exception: continue
    if name.startswith(('ctrl_', 'fib_', 'full_')):
        dec = S.get('decode', {})
        rec = {"exp": "B-reflsearch", "family": name, "d": None, "prime": S.get('p'), "result": f"genuine reflections: {dec.get('genuine')} (rational {dec.get('rational')}, quadratic {dec.get('quadratic')}, higher-degree solutions {dec.get('higher')}, spurious {dec.get('spurious')}); elim degree {S.get('degree')}", "rank": dec.get('genuine'), "seconds": S.get('seconds'), "job": "construct-rcr-chain1", "note": f"subst {S.get('subst')}; eigenvectors in {js}"}
    elif name.startswith('fermatcont'):
        rec = {"exp": "C-fermatcont", "family": f"base {S.get('base')}", "d": 6, "prime": S.get('p'), "result": f"orders: {[(o['order'], o['consistent']) for o in S.get('orders', [])]}; rank L {S.get('rankL')}, nullity {S.get('nullity')}", "rank": S.get('rankL'), "seconds": S.get('seconds'), "job": "construct-rcr-chain2", "note": "free variables set to 0 at each order"}
    elif name.startswith('gmrank'):
        rec = {"exp": "G-gm-rank", "family": "nabla-span of the (3,0)-forms in V_11", "d": None, "prime": S.get('p'), "result": f"ranks by derivative order: {S.get('ranks')} (basis {S.get('NB')})", "rank": (S.get('ranks') or [None])[-1], "seconds": S.get('seconds'), "job": "construct-rcr-chain1", "note": f"lambda0 {S.get('lam0')}, K {S.get('K')}, dims {S.get('dims')}"}
    elif name.startswith('critpts'):
        rec = {"exp": "F-critpts", "family": "critical points of Y0/Y1 on X", "d": None, "prime": S.get('p'), "result": f"elim degree {S.get('degree')} sqfree {S.get('sqfree')} factor degrees {S.get('factor_degrees')}; rational types {S.get('rational_types')}", "rank": S.get('degree'), "seconds": S.get('seconds'), "job": "construct-rcr-chain2", "note": f"n_rational {S.get('n_rational')}"}
    else:
        continue
    key = rec['exp'] + '|' + rec['family']
    if key in done: continue
    out.append(rec)
with open(L, 'a') as f:
    for rec in out: f.write(json.dumps(rec) + '\n')
print(f'appended {len(out)} records'); [print(json.dumps(r)[:300]) for r in out]
