"""(c): full E1+E2 sector-restricted test (with exact Farkas certificates) of non-Jacobian E1-kernel
representatives (rationally reconstructed, small integer coefficients) per G-sector at degree d.
usage: c_e1basis_farkas.py d [maxper] [random]"""
import sys, json, time, os, random
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import xlib
from xlib import *
F2, A, q = load_cycle_q2()
d = int(sys.argv[1]); maxper = int(sys.argv[2]) if len(sys.argv) > 2 else 999
mode = sys.argv[3] if len(sys.argv) > 3 else 'basis'
p = P2
JS = JacobianSectors(F2, p)
JS1 = JacobianSectors(F2, P1)
done = set()
for l in open(xlib.LOG):
    r = json.loads(l)
    if r.get('exp') == 'C-e1basis-farkas':
        done.add(r['D_hash'])
random.seed(d)
for key in all_sector_keys(JS.gens, d):
    reps, jd, kd = nonJ_representatives(JS, A, d, key)
    if not reps:
        continue
    # rational reconstruction; the kernel basis from nullspace is not canonical, so reconstruct each and
    # keep those that reconstruct with small coefficients; fall back to random small combinations verified mod P1
    cand = []
    for Dk in reps:
        Dq = poly_ratrec(Dk, p)
        if Dq is not None and max(abs(c) for c in Dq.values()) < 10 ** 6:
            cand.append(('rep', Dq))
    if mode == 'random' or not cand:
        # random small-integer combination inside the kernel: solve for small integer D directly over Q via nullspace lifting
        pass
    tested = 0
    for tag, Dq in cand[:maxper]:
        h = pstr_hash(Dq)
        if h in done:
            continue
        # verify E1 over a second prime
        kA = key_of_poly(A, JS1.gens)
        e1ok = express_in_J(JS1, pmod(pmul(A, Dq), P1), d + 12, key_add(key_of_poly(Dq, JS1.gens), kA)) is not None
        r = solve_fixed_D(F2, A, q, Dq)
        r.update(exp='C-e1basis-farkas', family=f'E1-nonJ-d{d}', sector=list(key), e1_modP1=e1ok, yparity=ypar(Dq), jdim=jd, kernel_dim=kd, n_nonJ=len(reps))
        if '_farkas_full' in r:
            fc = r['_farkas_full']
            json.dump(dict(D=to_str(Dq), q=q, F2=to_str(F2), A=to_str(A), H=int(r['H']), prime_used=r['primes'][-1], rows=fc['rows'], Y=fc['Y'], ell=fc['ell']),
                      open(f'/tmp/hodge/results/cycle-q2/wave2-excl/farkas/E1nonJ_d{d}_{h}.json', 'w'))
        print(d, key, to_str(Dq)[:80], r['status'], r['ranks'], r.get('farkas', {}).get('verified'), f"{r['seconds']}s", flush=True)
        xlib.log(r); tested += 1
print('done')
