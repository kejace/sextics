"""Logging helper for wave-2 excl lens (VM agent). Appends to wave2-excl.log.jsonl with engine 3 tag."""
import sys, time, json, os
sys.path.insert(0, '/tmp/hodge/results/dense-rational')
from dr import *
LOG2 = '/tmp/hodge/results/dense-rational/wave2-excl.log.jsonl'
OUT = '/tmp/hodge/results/dense-rational/wave2-excl'
os.makedirs(OUT, exist_ok=True)
def log2(**rec):
    rec.setdefault('lens', 'excl'); rec.setdefault('engine', 3); rec.setdefault('agent', 'vm'); rec.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    with open(LOG2, 'a') as f: f.write(json.dumps(rec) + '\n')
    return rec
def par_i(i): return 0 if i < 3 else 1
