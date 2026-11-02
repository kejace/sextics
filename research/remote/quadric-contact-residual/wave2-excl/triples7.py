"""Job excl-triples7: all 56 triples of (Z/2)^3 parity classes at degree 7 (~42-48 unknowns each), msolve -g 1 (threads from argv).
Usage: triples7.py threads [timeout]"""
import sys, itertools, subprocess, os
threads = sys.argv[1] if len(sys.argv) > 1 else '8'; tmo = sys.argv[2] if len(sys.argv) > 2 else '3600'
classes = [(b0, b2, y) for b0 in (0, 1) for b2 in (0, 1) for y in (0, 1)]
def conj(c): return f'b0={c[0]}&b2={c[1]}&y{c[2]}'
for tri in itertools.combinations(classes, 3):
    tag = 'd7_tri_' + '_'.join(''.join(map(str, c)) for c in tri)
    pred = '|'.join(conj(c) for c in tri)
    env = dict(os.environ, MSOLVE_TIMEOUT=tmo)
    r = subprocess.run(['python3', 'sub_closure.py', 'd7mixed.quads.local.json', tag, pred, threads], capture_output=True, text=True, env=env)
    for l in r.stderr.splitlines():
        if 'sub-closure-msolve' in l: print(l[:300], flush=True)
print('TRIPLES7 DONE')
