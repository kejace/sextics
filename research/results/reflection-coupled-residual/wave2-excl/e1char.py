"""E1-admissible denominators by Klein character: kernel of [A] on (S/J)_d^chi, with good-reduction checks.
Usage: python3 e1char.py <dmin> <dmax> <prime> [prime2]   (prime2: reconstruct kernel reps rationally when kernel != 0)
Receipts: wave2-excl/e1char.jsonl ; reps: wave2-excl/e1reps_d<d>_chi<ab>_p<prime>.json and e1reps_d<d>_chi<ab>_Q.json"""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
from hodge import core
dmin, dmax, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
p2 = int(sys.argv[4]) if len(sys.argv) > 4 else None
job = 'e1char'
for d in range(dmin, dmax + 1):
    for c in L.CHARS:
        tag = f'd{d}_chi{c[0]}{c[1]}'
        t0 = time.time()
        RT = L.CharReducer(d + 6, L.cadd(c, L.CHI_A), p)
        r = L.e1_kernel(d, c, p, red_target=RT)
        reps = r.pop('reps')
        json.dump([L.sym_str(P, p) for P in reps], open(f'{L.OUT}/e1reps_{tag}_p{p}.json', 'w'), indent=0)
        r['exp'] = 'e1-kernel-by-character'
        L.receipt(job, **r)
        if p2 and r['kernel_mod_J'] > 0:
            RT2 = L.CharReducer(d + 6, L.cadd(c, L.CHI_A), p2)
            r2 = L.e1_kernel(d, c, p2, red_target=RT2)
            reps2 = r2.pop('reps')
            r2['exp'] = 'e1-kernel-by-character-prime2'
            L.receipt(job, **r2)
            if len(reps2) == len(reps):
                rat = []
                allok = True
                for P1, P2 in zip(reps, reps2):
                    ok, Pq = L.ratrecon_polys([(p, P1), (p2, P2)])
                    allok &= ok
                    rat.append(core.to_str({e: v for e, v in Pq.items() if v is not None}) if ok else None)
                json.dump(rat, open(f'{L.OUT}/e1reps_{tag}_Q.json', 'w'), indent=0)
                L.receipt(job, exp='e1-kernel-ratrecon', d=d, chi=c, primes=[p, p2], n_classes=len(reps),
                          all_reconstructed=allok, note='canonical normal forms w.r.t. rref of J_d^chi; CRT of two primes')
            else:
                L.receipt(job, exp='e1-kernel-ratrecon', d=d, chi=c, primes=[p, p2], note='KERNEL DIM DISAGREES BETWEEN PRIMES',
                          dims=[len(reps), len(reps2)])
print('E1CHAR DONE', flush=True)
