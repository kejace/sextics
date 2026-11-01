import sys, json, time
from inv import *
d = int(sys.argv[1])
res, Ds, Ns, quads = quadratic_closure(d, verbose=True, out_prefix='/tmp/hodge/results/quadric-contact-residual/closure')
print(json.dumps(res))
