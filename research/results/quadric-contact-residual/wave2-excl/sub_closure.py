"""Sub-family closures of a quads JSON system: keep the unknowns whose character satisfies a predicate, export (c^2 syntax),
run msolve -g 1, and (optionally) the block Macaulay test at degree k.  Predicates: 'y0','y1' (tau_Y parity), 'b0=0','b0=1',
'b2=0','b2=1', 'h=k' (h1 character), combinations joined by '&'.  Usage: sub_closure.py quads.json tag pred threads [macaulay_k] [workers]"""
import sys, json, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msverdict import verdict
from macaulay import closure
fn, tag, pred = sys.argv[1], sys.argv[2], sys.argv[3]
threads = int(sys.argv[4]) if len(sys.argv) > 4 else 8
mk = int(sys.argv[5]) if len(sys.argv) > 5 else 0
workers = int(sys.argv[6]) if len(sys.argv) > 6 else 8
LOG = os.environ.get('HODGE_W2LOG', '/tmp/hodge/results/quadric-contact-residual/wave2-excl.log.jsonl')
def log(exp, **kw):
    rec = {'exp': exp, 'time': time.strftime('%Y-%m-%dT%H:%M:%S')}; rec.update(kw)
    open(LOG, 'a').write(json.dumps(rec) + '\n'); print('[LOG]', json.dumps(rec)[:500], file=sys.stderr, flush=True)
J = json.load(open(fn)); p = J['p']; chi = [tuple(c) for c in J['chi']]
def ok1(c, conj):
    for t in conj.split('&'):
        if t == 'y0' and c[3] != 0: return False
        if t == 'y1' and c[3] != 1: return False
        if t.startswith('b0=') and c[0] != int(t[3:]): return False
        if t.startswith('b2=') and c[1] != int(t[3:]): return False
        if t.startswith('h=') and c[2] != int(t[2:]): return False
    return True
def ok(c):
    return any(ok1(c, conj) for conj in pred.split('|'))
keep = [a for a in range(len(chi)) if ok(chi[a])]
quads = [{(a, b): v for a, b, v in cf} for cf in J['quads']]
ks = set(keep)
Q = []
for cf in quads:
    cf2 = {(a, b): v for (a, b), v in cf.items() if a in ks and b in ks}
    if cf2: Q.append(cf2)
prefix = os.path.join(os.path.dirname(os.path.abspath(fn)), f'sub_{tag}')
with open(prefix + '.ms', 'w') as f:
    eqs = ['+'.join((f'{v}*c{a}^2' if a == b else f'{v}*c{a}*c{b}') for (a, b), v in cf.items()) for cf in Q]
    f.write(','.join(f'c{a}' for a in keep) + '\n' + str(p) + '\n' + ',\n'.join(eqs) + '\n')
json.dump(keep, open(prefix + '.subset.json', 'w'))
t0 = time.time()
info = {'tag': tag, 'pred': pred, 'n_unknowns': len(keep), 'n_quadrics': len(Q), 'source': os.path.basename(fn)}
if threads > 0:
    v = verdict(prefix + '.ms', threads=threads, timeout=int(os.environ.get('MSOLVE_TIMEOUT', '14400')), verbose=1, logfile=prefix + '.msolve.log')
    info['msolve'] = {k: val for k, val in v.items() if k not in ('log_tail', 'pure_powers')}
    info['result'] = 'EXCLUDED (only c=0)' if v.get('verdict', '').startswith('ONLY-ORIGIN') else v.get('verdict')
    log('sub-closure-msolve', seconds=round(time.time() - t0, 1), **info)
if mk:
    t1 = time.time()
    s = closure(fn, mk, workers, keep, log=lambda m: print(m, file=sys.stderr, flush=True))
    json.dump(s, open(prefix + f'.macaulay{mk}.json', 'w'))
    log('sub-closure-macaulay', tag=tag, pred=pred, seconds=round(time.time() - t1, 1), **{k: v for k, v in s.items() if k != 'blocks'},
        result='EXCLUDED (I_k = S_k)' if s['I_k_equals_S_k'] else 'not closed at this degree')
print('DONE', file=sys.stderr)
