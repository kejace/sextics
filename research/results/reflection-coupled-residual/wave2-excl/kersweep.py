"""Fixed-D tests (rigorfix protocol) of E1-kernel representatives and random combinations.
Usage: python3 kersweep.py <d> <chi as ab> [nrandom] [--full]
Reads wave2-excl/e1reps_d<d>_chi<ab>_Q.json (rational reps; falls back to the mod-p file)."""
import sys, json, random, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
import rigorfix
from hodge import core
from fractions import Fraction
d = int(sys.argv[1]); ab = sys.argv[2]
nrand = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith('--') else 5
full = '--full' in sys.argv
fn = f'{L.OUT}/e1reps_d{d}_chi{ab}_Q.json'
try:
    reps = json.load(open(fn))
    src = 'rational'
except Exception:
    fn = f'{L.OUT}/e1reps_d{d}_chi{ab}_p2147483587.json'
    reps = json.load(open(fn)); src = 'modp(2147483587)'
reps = [r for r in reps if r]
print(f'{len(reps)} reps from {fn} ({src})', flush=True)
random.seed(d * 100 + int(ab))
for k, r in enumerate(reps):
    try:
        rigorfix.run(f'e1ker-d{d}-chi{ab}-basis{k}', r, full=full)
    except Exception as ex:
        L.receipt('rigorfix', exp='rigorfix', family=f'e1ker-d{d}-chi{ab}-basis{k}', error=repr(ex))
polys = [core.parse(r) for r in reps]
for t in range(nrand):
    coeffs = [random.randint(-3, 3) for _ in polys]
    D = {}
    for c, P in zip(coeffs, polys):
        if c: D = core.add(D, core.scale(P, Fraction(c)))
    if not D: continue
    expr = core.to_str(D)
    try:
        rigorfix.run(f'e1ker-d{d}-chi{ab}-rand{t}', expr, full=False)
    except Exception as ex:
        L.receipt('rigorfix', exp='rigorfix', family=f'e1ker-d{d}-chi{ab}-rand{t}', error=repr(ex))
print('KERSWEEP DONE', flush=True)
