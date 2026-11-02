#!/usr/bin/env python3
"""Calibration: success rate of the LM VarPro search on the Fermat sextic (known solutions f6*L) at d = 6,7,8, (0,0)."""
import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import json, time, math
import multiprocessing as mp
import numpy as np
import w2lib
from hodge import core
F, A, q = w2lib.load()
F0 = core.parse('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
for d in (6, 7, 8):
    basis, gauge, info = w2lib.build_family(F0, A, q, d, (0, 0), 'colon')
    fam = w2lib.SectorVarPro(F0, A, q, basis, (0, 0), gauge=gauge)
    t0 = time.time()
    def run(seed):
        r = fam.lm_from(seed, maxiter=100)
        return (seed, r['relres'], r['iters'], r['status'], r['time'])
    with mp.get_context('fork').Pool(workers) as pool:
        res = pool.map(run, range(n))
    rels = sorted(x[1] for x in res)
    succ = sum(1 for x in rels if x < 1e-10)
    rec = {'ts': w2lib.now(), 'step': 'control-fermat', 'd': d, 'sector': [0, 0], 'K': fam.K, 'gauge': fam.ng, 'dimP': fam.dimP,
           'rows': fam.R, 'nN': fam.nN, 'nQ': fam.nQ, 'restarts': n, 'successes': succ, 'best': rels[0], 'median': rels[len(rels)//2],
           'mean_time': round(float(np.mean([x[4] for x in res])), 2), 'wall': round(time.time() - t0, 1),
           'minima': sorted({f'{x:.3e}' for x in rels})[:10]}
    print(json.dumps(rec), flush=True)
    w2lib.append_json('wave2-search.remote.log.jsonl', rec)
