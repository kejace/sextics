"""Test the invariant E1-passing space K_d^G = (J:A)_d^G with the invariant solver (mod p consistency),
on kernel representatives, random points of K, and structured guesses."""
import sys, time, json, pickle, random
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
from h import *
from invsolve import *
F, A, q = load()
p = P0
mul, add, sub, sc, pw = core.mul, core.add, core.sub, core.scale, core.power
Fi = [core.pderiv(F, i) for i in range(6)]

def J_inv_basis(d, pD):
    out = []
    for slot in (0, 3):
        par_m = (pD - (0 if slot == 0 else 1)) % 2
        for m in orbit_reps(d - 5, par_m, lambda E: canon_H(E, slot)):
            out.append(sym(core.shift(Fi[slot], m, 1)))
    return out

def modp_lift(Pol):
    return {e: (v % p if v % p < p // 2 else v % p - p) for e, v in core.reduce_mod(Pol, p).items()}

def run(D, lab, fam, note=''):
    t0 = time.time()
    if not is_invariant(D):
        print(lab, 'not invariant'); return None
    ok, info = e1_inv(F, A, q, D)
    if not ok:
        log(lab, core.degree(D), fam, False, None, info['rank'], info['cols'] - info['rank'], time.time() - t0, note + f' | inv E1 fail {info} engine=inv-v1')
        print(f'{lab:40s} E1=False', flush=True); return False
    S = InvSystem(F, A, q, D, p=p)
    r = S.solve()
    log(lab, core.degree(D), fam, True, r['consistent'], r['rank'], r['nullity'], time.time() - t0,
        note + f" | inv full modp rows={r['rows']} cols={r['cols']} engine=inv-v1")
    print(f"{lab:40s} d={core.degree(D)} E1=True consistent={r['consistent']} rank={r['rank']} null={r['nullity']} rows={r['rows']} cols={r['cols']} ({time.time()-t0:.1f}s)", flush=True)
    if r['consistent']:
        print('!!!!! CONSISTENT mod p:', lab, flush=True)
        pickle.dump(D, open(f'{RES}/consistent_{lab.replace(" ","_").replace("/","_")}.pkl', 'wb'))
    return r['consistent']

if __name__ == '__main__':
    random.seed(1)
    degs = [(8, 0), (9, 0), (9, 1)] if len(sys.argv) < 2 else [tuple(map(int, a.split(','))) for a in sys.argv[1:]]
    for d, pD in degs:
        info, polys = pickle.load(open(f'{RES}/e1inv_d{d}_p{pD}.pkl', 'rb'))
        Jb = J_inv_basis(d, pD)
        print(f'=== degree {d} parity {pD}: kernel reps {len(polys)}, J_d^G dim {len(Jb)}', flush=True)
        for i, Pk in enumerate(polys):
            run(Pk, f'K{d}p{pD}_rep{i}', f'E1-kernel-inv-d{d}', note='kernel representative mod p as returned by rref')
        # random points in the quotient-span (kernel reps only)
        for t in range(4):
            D = {}
            for Pk in polys:
                D = add(D, sc(Pk, random.randrange(1, p)))
            run(modp_lift(D), f'K{d}p{pD}_randker{t}', f'E1-kernel-inv-d{d}', note='random combination of kernel reps (mod p)')
        # random points in the full K = ker + J
        for t in range(6):
            D = {}
            for Pk in polys:
                D = add(D, sc(Pk, random.randrange(1, p)))
            for Jk in Jb:
                D = add(D, sc(Jk, random.randrange(0, p)))
            run(modp_lift(D), f'K{d}p{pD}_randfull{t}', f'E1-kernel-inv-d{d}', note='random point of (J:A)_d^G incl. J part (mod p)')
        # pure J elements (random) as control
        for t in range(3):
            D = {}
            for Jk in Jb:
                D = add(D, sc(Jk, random.randrange(0, p)))
            run(modp_lift(D), f'J{d}p{pD}_rand{t}', f'J-inv-d{d}', note='random invariant element of J_d (mod p)')
