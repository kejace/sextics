import sys, json, time
from inv import *
d = int(sys.argv[1])
if len(sys.argv) > 2:
    c = tuple(int(x) for x in sys.argv[2].split(','))
    res, Ds, Ns, quads = quadratic_closure(d, verbose=True, out_prefix='/tmp/hodge/results/quadric-contact-residual/closure', c=c, sym=False)
else:
    res, Ds, Ns, quads = quadratic_closure(d, verbose=True, out_prefix='/tmp/hodge/results/quadric-contact-residual/closure')
print(json.dumps(res))
