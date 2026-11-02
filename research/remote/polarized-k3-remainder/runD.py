"""Run the fixed-D solver on one denominator given as an expression in the (U,V) frame or (X,Y) frame.
usage: runD.py '<family>' '<label>' '<expr>' [uv|xy] [fast|dense] [target]
In the uv frame, the expression is written in variables X0..X2 (= U) and Y0..Y2 (= V) and is converted
to the physical X,Y frame by U = X+Y, V = X-Y."""
import sys, time, json, resource
sys.path.insert(0, '/home/kejace/devel/sextics-work'); sys.path.insert(0, '/tmp/hodge/results/polarized-k3-remainder')
from hodge import core, fixed_d
from logutil import log
from fractions import Fraction as Fr

def uv_to_xy(P):
    def var(i):
        e = [0]*6; e[i] = 1; return {tuple(e): Fr(1)}
    X = [var(i) for i in range(3)]; Y = [var(3+i) for i in range(3)]
    imgs = [core.add(X[i], Y[i]) for i in range(3)] + [core.sub(X[i], Y[i]) for i in range(3)]
    out = {}
    for e, c in P.items():
        term = {(0,)*6: Fr(c)}
        for i in range(6):
            if e[i]:
                term = core.mul(term, core.power(imgs[i], e[i]))
        out = core.add(out, term)
    return out

def main():
    family, label, expr = sys.argv[1], sys.argv[2], sys.argv[3]
    frame = sys.argv[4] if len(sys.argv) > 4 else 'uv'
    mode = sys.argv[5] if len(sys.argv) > 5 else 'fast'
    want_target = len(sys.argv) > 6 and sys.argv[6] == 'target'
    F, A, q, rec = core.load_instance('polarized-k3-remainder')
    D = core.parse(expr)
    if frame == 'uv':
        D = uv_to_xy(D)
    d = core.degree(D)
    par = core.parity_of(D)
    print(f'[{label}] d={d} terms={len(D)} parity={par}', flush=True)
    if par is None:
        log(label, d, family, None, None, -1, -1, 0.0, 'mixed parity in XY frame; skipped')
        return
    t0 = time.time()
    import os
    wl = open('/tmp/hodge/results/polarized-k3-remainder/target_whitelist.txt').read().split() if os.path.exists('/tmp/hodge/results/polarized-k3-remainder/target_whitelist.txt') else []
    want_target = want_target and (label in wl)
    if d >= 8 and os.environ.get('FORCE_E1') is None:
        ok, info = True, {'skipped': True}
        print(f'[{label}] E1 test skipped at d>=8 (full system decides)', flush=True)
    else:
        ok, info = fixed_d.e1_solvable(F, A, q, D)
    t1 = time.time()
    print(f'[{label}] E1 solvable={ok} {info} ({t1-t0:.1f}s)', flush=True)
    if not ok:
        log(label, d, family, False, False, info['rank'], info['cols'] - info['rank'], t1 - t0,
            f'E1 fails (rank {info["rank"]}, aug {info["rank_aug"]}); frame={frame}; expr={expr}')
        return
    if mode == 'fast':
        out = fixed_d.solve_fixed_D_fast(F, A, q, D, verbose=True, want_exact=True)
        t2 = time.time()
        cons = out['consistent']
        note = f'fast; rows={out["rows"]} cols={out["cols"]}; frame={frame}; expr={expr}'
        if cons:
            note += f'; exact={out.get("exact")}'
            if out.get('exact'):
                w = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}
                json.dump(w, open(f'/tmp/hodge/results/polarized-k3-remainder/witness-candidate-{label}.json', 'w'))
                print('WITNESS FOUND', flush=True)
        log(label, d, family, True, cons, out['rank'], out['nullity'], t2 - t0, note)
        print(f'[{label}] consistent={cons} rank={out["rank"]} nullity={out["nullity"]} ({t2-t0:.1f}s) maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024} MB', flush=True)
    else:
        status, data = fixed_d.solve_fixed_D(F, A, q, D, verbose=True)
        t2 = time.time()
        cons = status == 'solved'
        note = f'dense status={status}; frame={frame}; expr={expr}'
        if cons:
            w = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in data['N']], 'Q': core.to_str(data['Q'])}
            json.dump(w, open(f'/tmp/hodge/results/polarized-k3-remainder/witness-candidate-{label}.json', 'w'))
            print('WITNESS FOUND', flush=True)
        log(label, d, family, True, cons, data.get('rank'), data.get('nullity'), t2 - t0, note)
        print(f'[{label}] status={status} {dict((k,v) for k,v in data.items() if k not in ("N","Q"))} ({t2-t0:.1f}s)', flush=True)
    if want_target:
        t3 = time.time()
        dim, basis, S = fixed_d.target_image(F, q, D, verbose=True)
        json.dump({'label': label, 'D_xy': core.to_str(D), 'dim': dim, 'basis': [core.to_str(b) for b in basis]},
                  open(f'/tmp/hodge/results/polarized-k3-remainder/target-image-{label}.json', 'w'))
        print(f'[{label}] target_image dim={dim} ({time.time()-t3:.1f}s)', flush=True)

main()
