"""Re-verification of the E1 (annihilator) profile at a new prime with independent code (w2.e1_family).
For each d, parity: dim K_d = {D : A D in J}, dim J_d, kernel mod J; plus Fermat-value check of rank J_{d+6}."""
import sys, json, time
sys.path.insert(0, '.')
import w2
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000007
ds = [int(x) for x in sys.argv[2:]] or list(range(1, 10))
R, Fp, Ap, Fi, F, A = w2.load(p)

def fermat_hf(k, par):
    # standard monomials of the Fermat ring: exponents 0..4, degree k, parity par
    import itertools
    cnt = 0
    for e in w2.monomials(k, par):
        if max(e) <= 4:
            cnt += 1
    return cnt

for d in ds:
    for pD in (0, 1):
        t0 = time.time()
        fam = w2.e1_family(R, Fp, Ap, Fi, d, pD, verbose=False)
        rJ, _ = w2.jacobian_dim(R, Fi, d, pD)
        prow = (1 + pD) % 2
        rJ6, _ = w2.jacobian_dim(R, Fi, d + 6, prow)
        nS = len(w2.monomials(d + 6, prow))
        fermat = nS - fermat_hf(d + 6, prow)
        nSd = len(w2.monomials(d, pD)); fermat_d = nSd - fermat_hf(d, pD)
        rec = dict(exp='e1-kernel', p=p, d=d, parity=pD, dimK=fam['dimK'], dimJ_d=rJ, kernel_mod_J=fam['dimK'] - rJ,
                   rank_J_d_plus_6=rJ6, fermat_rank_J_d_plus_6=fermat, rank_matches_fermat=(rJ6 == fermat),
                   rank_J_d_matches_fermat=(rJ == fermat_d), rows=fam['rows'], cols=fam['cols'],
                   seconds=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True)
        w2.log(rec)
        if fam['dimK'] - rJ > 0:
            json.dump({'p': p, 'd': d, 'parity': pD, 'basis': [[str(D), [str(x) for x in n]] for D, n in fam['basis']]},
                      open(f'e1fam_d{d}_p{pD}_{p}.json', 'w'))
