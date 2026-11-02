"""Rigorous E1 screen of structured candidates: A*D in J(F) mod p?  With smooth reduction mod p, 'not in J mod p'
implies 'not in J over Q' (rank argument), hence D (and every power/multiple of it in the Jacobian-ring sense) is
excluded.  usage: python3 e1_screen.py candidates.jsonl [p]"""
import sys, json, time
from common import *
from hodge import jacobian
from hodge.core import mul, reduce_mod
F, A, q, rec = core.load_instance(IID)
p = int(sys.argv[2]) if len(sys.argv) > 2 else core.PRIMES[0]
for line in open(sys.argv[1]):
    c = json.loads(line); D = core.parse(c['D']); t0 = time.time()
    AD = reduce_mod(mul(A, D), p)
    ok, info = jacobian.in_jacobian(reduce_mod(F, p), AD, p)
    log({'exp': 'e1-screen', 'name': c['name'], 'd': c['d'], 'parityD': c['parity'], 'p': p, 'e1_in_J_modp': bool(ok),
         'rigorous_exclusion_over_C': (not ok), 'seconds': round(time.time() - t0, 1)})
