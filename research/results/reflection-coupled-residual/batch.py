"""Batch runner: for each (family, D-expression) compute served space, log, accumulate even-sector span."""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual')
from rcr import *
import fcntl

def run_batch(cands, span_file='span_even.json', tag=''):
    J = get_J6()
    av = A_vec()
    for fam, expr in cands:
        D = core.parse(expr)
        d = core.degree(D)
        if core.parity_of(D) is None:
            print('skip mixed parity', expr); continue
        t0 = time.time()
        try:
            ok1, info1 = fixed_d.e1_solvable(F, A, q, D, p=P)
        except Exception as ex:
            ok1, info1 = None, str(ex)
        res = served(D)
        # does D serve A directly?
        direct = res['span'].contains(av)
        # accumulate even-sector span (shared file, locked)
        with open('span.lock', 'w') as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            sp = Span.load(span_file) if os.path.exists(span_file) else Span(J.quot_dim)
            before = sp.dim()
            for piv, r in res['span_even'].rows:
                sp.add(r)
            after = sp.dim()
            inA = sp.contains(av)
            sp.save(span_file)
            fcntl.flock(lk, fcntl.LOCK_UN)
        secs = time.time() - t0
        note = (f"served dim_modJ={res['dim_modJ']} even={res['dim_even']} odd={res['dim_odd']} raw={res['dim_raw']}; "
                f"cum even span {before}->{after}; A in cum span: {inA}; rows={res['rows']} cols={res['cols']}")
        print(f"[{tag}] {fam:28s} d={d} e1={ok1} direct={direct} | {note} ({secs:.1f}s)", flush=True)
        log(expr, d, fam, ok1, direct, res['rows'] and None, None, secs, note,
            served_modJ=res['dim_modJ'], served_even=res['dim_even'], served_odd=res['dim_odd'], cum_even=after, A_in_cum=inA)
        if direct:
            print('!!!!!!!! DIRECT HIT: D =', expr, flush=True)
            with open('HITS.txt', 'a') as f: f.write(expr + '\n')
        if inA:
            print('!!!!!!!! A IN CUMULATIVE SPAN after', expr, flush=True)
            with open('HITS.txt', 'a') as f: f.write('CUMULATIVE: ' + expr + '\n')

if __name__ == '__main__':
    cands = [tuple(l.rstrip('\n').split('\t')) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
    run_batch(cands, tag=sys.argv[2] if len(sys.argv) > 2 else os.path.basename(sys.argv[1]))
