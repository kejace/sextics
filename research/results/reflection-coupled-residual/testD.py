"""Full fixed-D test (engine v2 fast compressed solver) + log. Usage: python3 testD.py <family> <D-expr> [--dense]"""
import sys, time, json
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
IID = 'reflection-coupled-residual'; LOG = f'/tmp/hodge/results/{IID}/log.jsonl'
F, A, q, rec = core.load_instance(IID)
fam, expr = sys.argv[1], sys.argv[2]
dense = '--dense' in sys.argv
D = core.parse(expr); d = core.degree(D)
t0 = time.time()
ok1, info1 = fixed_d.e1_solvable(F, A, q, D)
t1 = time.time()
print(f'[{fam}] d={d} e1={ok1} ({t1-t0:.1f}s) info={str(info1)[:200]}', flush=True)
ent = {"D": expr, "d": d, "family": fam, "e1": bool(ok1), "consistent": None, "rank": None, "nullity": None, "seconds": 0, "note": "", "engine": 2}
if ok1:
    try:
        if dense:
            st, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True)
            cons = st != 'inconsistent'
            ent.update(consistent=cons, rank=data.get('rank'), nullity=data.get('nullity'), note=f'dense status={st}')
            if st == 'solved':
                N, Q = data['N'], data['Q']
        else:
            out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True)
            cons = out['consistent']
            ent.update(consistent=bool(cons), rank=out.get('rank'), nullity=out.get('nullity'),
                       note=f"fast rows={out.get('rows')} cols={out.get('cols')} exact={out.get('exact')}")
            if cons and out.get('exact'):
                N, Q = out['N'], out['Q']
        if ent['consistent'] and 'N' in dir():
            okv, infov = core.verify_witness(F, A, q, D, N, Q)
            ent['note'] += f' VERIFIED={okv}'
            if okv:
                json.dump({'D': expr, 'N': [core.to_str(n) for n in N], 'Q': core.to_str(Q)},
                          open(f'/tmp/hodge/results/{IID}/witness-candidate.json', 'w'), indent=1)
                print('!!!!!!!! WITNESS FOUND AND VERIFIED', flush=True)
    except Exception as ex:
        ent['note'] = f'ERROR {ex!r}'
        print('ERROR', repr(ex), flush=True)
else:
    ent['note'] = 'E1 fails (A*D not in J); full system skipped'
ent['seconds'] = round(time.time() - t0, 1)
print(json.dumps(ent)[:600], flush=True)
with open(LOG, 'a') as f: f.write(json.dumps(ent) + '\n')
