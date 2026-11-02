"""Shared helpers for the wave-2 exclusion lens (dense-rational). Not part of the repo."""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *          # core, fixed_d, NV, monomials, shift, mul, add, sub, scale, pderiv, power, nmod_mat, load, hilbert_counts
LOG2 = '/tmp/hodge/results/dense-rational/wave2-excl.log.jsonl'
OUT = '/tmp/hodge/results/dense-rational/wave2-excl'
os.makedirs(OUT, exist_ok=True)
def log2(**rec):
    rec.setdefault('lens', 'excl'); rec.setdefault('engine', 2); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG2, 'a') as f: f.write(json.dumps(rec) + '\n')
    return rec
def par_i(i): return 0 if i < 3 else 1
