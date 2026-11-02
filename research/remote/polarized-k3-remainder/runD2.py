"""usage: runD2.py '<family>' '<label>' '<expr>' [uv|xy]   -- lean fixed-D solve, logs to log.jsonl.
E1 is skipped for d>=7 (the full system decides and E1 alone costs 10 min at d=8)."""
import sys, time, json, resource, os
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
            if e[i]: term = core.mul(term, core.power(imgs[i], e[i]))
        out = core.add(out, term)
    return out
from fastsolve import solve_lean
def main():
    family, label, expr = sys.argv[1], sys.argv[2], sys.argv[3]
    frame = sys.argv[4] if len(sys.argv) > 4 else 'uv'
    F, A, q, rec = core.load_instance('polarized-k3-remainder')
    D = core.parse(expr)
    if frame == 'uv': D = uv_to_xy(D)
    d = core.degree(D); par = core.parity_of(D)
    print(f'[{label}] d={d} terms={len(D)} parity={par}', flush=True)
    if par is None:
        log(label, d, family, None, None, -1, -1, 0.0, f'mixed Y-parity in XY frame; skipped; expr={expr}'); sys.exit()
    # legality guard: D must not be divisible by F (mod-p rank test, exact enough to flag)
    if d >= 6:
        from flint import nmod_mat
        pp = core.PRIMES[0]; Fp = core.reduce_mod(F, pp); Dp = core.reduce_mod(D, pp)
        gm = core.monomials(d - 6); colsF = [core.shift(Fp, m, 1, pp) for m in gm]
        rows = {}
        for c in colsF + [Dp]:
            for e in c: rows.setdefault(e, len(rows))
        R, Cn = len(rows), len(colsF)
        flat = [0] * (R * (Cn + 1))
        for j, c in enumerate(colsF):
            for e, v in c.items(): flat[rows[e] * (Cn + 1) + j] = v
        for e, v in Dp.items(): flat[rows[e] * (Cn + 1) + Cn] = v
        Maug = nmod_mat(R, Cn + 1, flat, pp); M0 = nmod_mat(R, Cn, [flat[r * (Cn + 1) + j] for r in range(R) for j in range(Cn)], pp)
        if Maug.rank() == M0.rank():
            print(f'[{label}] ILLEGAL: D is divisible by F (mod p)', flush=True)
            log(label, d, family, None, None, -1, -1, 0.0, f'ILLEGAL: F | D; skipped; expr={expr}'); sys.exit()
    t0 = time.time()
    e1 = None
    if d <= 6:
        ok, info = fixed_d.e1_solvable(F, A, q, D); e1 = ok
        print(f'[{label}] E1 solvable={ok} {info} ({time.time()-t0:.1f}s)', flush=True)
        if not ok:
            log(label, d, family, False, False, info['rank'], info['cols']-info['rank'], time.time()-t0, f'E1 fails; frame={frame}; expr={expr}'); sys.exit()
    out = solve_lean(F, A, q, D, verbose=True)
    t2 = time.time()
    cons = out['consistent']
    note = f'lean; rows={out["rows"]} cols={out["cols"]} k={out["k"]} t_build={out["t_build"]} t_compress={out["t_compress"]} t_rref={out["t_rref"]}; frame={frame}; expr={expr}'
    if cons:
        note += f'; exact={out.get("exact")}'
        if out.get('exact'):
            w = {'d': d, 'D': core.to_str(D), 'N': [core.to_str(n) for n in out['N']], 'Q': core.to_str(out['Q'])}
            okv, infov = core.verify_witness(F, A, q, D, out['N'], out['Q'])
            w['verify'] = [okv, str(infov)]
            json.dump(w, open(f'/tmp/hodge/results/polarized-k3-remainder/witness-candidate-{label}.json', 'w'))
            print('WITNESS FOUND verify=', okv, flush=True)
    log(label, d, family, True if e1 is None else e1, cons, out['rank'], out['nullity'], t2-t0, note)
    print(f'[{label}] consistent={cons} rank={out["rank"]} nullity={out["nullity"]} ({t2-t0:.1f}s) maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024} MB', flush=True)

if __name__ == '__main__':
    main()
