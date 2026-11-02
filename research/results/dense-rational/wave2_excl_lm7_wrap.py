#!/usr/bin/env python3
"""Run wave-1 lm_mixed7.py (degree-7 Hessian relaxation, float LM) and copy its JSON summary into the wave-2 log.
Usage: python3 wave2_excl_lm7_wrap.py nstarts seed sector"""
import sys, subprocess, json
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from wave2_excl_common import *
nstarts, seed, sector = sys.argv[1], sys.argv[2], sys.argv[3]
r = subprocess.run(['python3', 'lm_mixed7.py', nstarts, seed, sector], capture_output=True, text=True, cwd='/tmp/hodge/results/dense-rational')
last = [l for l in r.stdout.splitlines() if l.startswith('{')]
rec = json.loads(last[-1]) if last else {'error': r.stderr[-500:]}
rec.update({'exp': 'lm-hessian-relaxation-d7', 'd': 7, 'script': 'lm_mixed7.py (wave 1)'})
if 'best10' in rec: rec['result'] = 'no near-solution (all scale-free residuals >= 1e-8)' if rec['best10'][0] >= 1e-8 else 'NEAR-SOLUTION FOUND'
log2(**rec); print(json.dumps(rec))
