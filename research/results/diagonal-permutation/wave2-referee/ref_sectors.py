"""(tau x (01))-sector dimensions of the Jacobian ring from the equivariant Koszul series, and forced E1 kernels.
[A]. maps sector (par, sig) of R_d to (1-par, -sig) of R_{d+6}  (A is Y-odd and (01)-anti-invariant)."""
import sys, json; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
T = R.sector_series(30)
out = {'table': {}, 'forced_kernels': {}}
hdr = 'n   ' + ' '.join(f'{n:>5}' for n in range(6, 19)); print(hdr)
for par, sig, lab in ((0, 1, '+,+'), (0, -1, '+,-'), (1, 1, '-,+'), (1, -1, '-,-')):
    row = [T[(n, par, sig)] for n in range(6, 19)]
    out['table'][lab] = row
    print(f'{lab} ' + ' '.join(f'{v:>5}' for v in row))
for d in range(1, 13):
    fk = {}
    for par in (0, 1):
        for sig in (1, -1):
            src = T[(d, par, sig)]; tgt = T[(d + 6, 1 - par, -sig)]
            fk[f'par{par}_sig{sig:+d}'] = {'src': src, 'tgt': tgt, 'forced_min_kernel': max(0, src - tgt)}
    out['forced_kernels'][d] = fk
    print('d', d, {k: v['forced_min_kernel'] for k, v in fk.items()})
# cross-check: full-parity dims equal sums of sectors
plus, minus = R.koszul_series(30)
assert all(T[(n, 0, 1)] + T[(n, 0, -1)] == plus[n] and T[(n, 1, 1)] + T[(n, 1, -1)] == minus[n] for n in range(25))
print('sector sums match parity-refined Koszul series: OK')
json.dump(out, open('ref_sectors.json', 'w'), indent=1)
