"""Calibration: on the split sextic F0 = f + g (t=0), D = f should serve the whole bidegree-(3,3) sector (100 classes + J_6)."""
import time, json
import lib
from lib import *
F0 = core.add(f, g)
t0 = time.time()
dim, basis, S = fixed_d.target_image(F0, 1, f)
print('t=0: served dim by D=f on F0=f+g:', dim, f'({time.time()-t0:.0f}s)')
# how many of the served classes are of block bidegree (3,3)?  count basis elements' bidegree support
bideg = {}
for b in basis:
    key = tuple(sorted({(sum(e[:3]), sum(e[3:])) for e in b}))
    bideg[key] = bideg.get(key, 0) + 1
print('bidegree supports of served basis:', bideg)
ok, info = fixed_d.e1_solvable(F0, A, 1, f)
print('E1 for A=X0^3Y0^3 over D=f at t=0:', ok)
status, data = fixed_d.solve_fixed_D(F0, A, 1, f)
print('fixed-D solve at t=0 with D=f:', status, {k: v for k, v in data.items() if k not in ('N', 'Q')})
if status == 'solved':
    print('N =', [core.to_str(n) for n in data['N']]); print('Q =', core.to_str(data['Q']))
log({'D': 'f on F0=f+g (t=0 calibration)', 'd': 6, 'family': 'calibration-t0', 'e1': bool(ok), 'consistent': status == 'solved',
     'rank': dim, 'nullity': dim - 18, 'seconds': round(time.time() - t0, 1), 'note': f'served dim {dim} on split F0; status {status}', 'engine': 2})
