"""Full fixed-D solves on the remote box for a list of D's stored as a pickle (list of (name, family, note, poly-mod-p or exact poly))."""
import sys, time, json, pickle
sys.path.insert(0,'/tmp/hodge/results/eight-plane-configuration')
from lib import *
import lib
_log = lib.log
def log2(e): e['engine'] = 2; e['host'] = 'remote'; _log(e)
lib.log = log2
items = pickle.load(open(sys.argv[1], 'rb'))
p_fam = int(sys.argv[2]) if len(sys.argv) > 2 else None   # prime the D-coefficients live in (None = exact rationals)
for name, fam, note, D in items:
    if p_fam:
        # lift symmetric representatives to integers (coefficients mod p -> centred)
        D = {e: Fraction(v if v <= p_fam//2 else v - p_fam) for e, v in D.items()}
    e = test_D(D, name, fam, note=note + (f' (coeffs taken mod {p_fam}, centred lift; solve mod engine primes)' if p_fam else ''), fast=True, want_exact=True, p=p_fam)
    print(name, {k: e.get(k) for k in ('d','e1','consistent','rank','nullity','rows','cols','seconds','note')}, flush=True)
