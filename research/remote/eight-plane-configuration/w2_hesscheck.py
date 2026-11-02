import sys; sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from w2lib import *
for pD in (0, 1):
    res = hess_quadrics(6, pD, p=PRIME, verbose=True)
    rows = quadric_rows(res)
    nW = res['nW']; Wm = res['Wm']
    # Euler vector: w_i = z_i  -> coordinate 1 on (i, e_i)
    ev = [0] * nW
    for k, (i, m) in enumerate(Wm):
        if m[i] == 1 and sum(m) == 1:
            ev[k] = 1
    vals = [sum(v * ev[k] * ev[l] for (k, l), v in row.items()) % PRIME for row in rows]
    print('pD', pD, 'rank', res['rank_quadrics'], 'Euler in family:', sum(ev), 'nonzero quadrics at Euler:', sum(1 for v in vals if v))
    write_msolve(res, f'{WD}/hess6fix-quadrics-d6-p{pD}.ms')
