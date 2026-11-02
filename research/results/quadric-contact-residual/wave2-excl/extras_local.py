"""Local sample of the extras screen: sector-restricted full fixed-D solves of the non-Jacobian E1-kernel basis elements of one sector."""
import sys, time, json
from w2lib import *
from family import e1_profile_sector
d = int(sys.argv[1]); c = tuple(int(x) for x in sys.argv[2]); nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 100
prof, extras = e1_profile_sector(d, c)
print('sector', c, prof, 'extras', len(extras), file=sys.stderr)
ncons = 0
for i, P in enumerate(extras[:nmax]):
    t0 = time.time()
    S = SectorSystem(P, verbose=False); res = S.solve()
    ok = None
    if res['consistent']:
        N, Q = S.unpack(res['x']); E = core.defect(Fp, Ap, P, 1, N, p); ok = (sub(E, mul(power(Fp, 2, p), Q, p), p) == {}); ncons += 1
        json.dump({'D': core.to_str(P), 'N': [core.to_str(x) for x in N], 'Q': core.to_str(Q), 'p': p}, open(PK + f'/CONSISTENT_d{d}_{"".join(map(str,c))}_{i}.json', 'w'))
    log('extras-screen', d=d, chi=str(c), i=i, consistent=res['consistent'], verified_mod_p=ok, rank=res['rank'], nullity=res['nullity'], rows=res['rows'], cols=res['cols'],
        seconds=round(time.time() - t0, 1), result='CONSISTENT' if res['consistent'] else 'inconsistent', note='local sample; basis element of the non-Jacobian E1 kernel (screen, not a family exclusion)')
    print(i, res['consistent'], res['rank'], res['nullity'], f'{time.time()-t0:.1f}s', flush=True)
print('consistent:', ncons, 'of', min(len(extras), nmax))
