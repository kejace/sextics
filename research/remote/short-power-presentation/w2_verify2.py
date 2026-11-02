import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
F, A, q, rec = load()
D7 = core.parse(json.load(open('kernel_d7_p0.json'))['basis_rational'][0])
D7o = core.parse(json.load(open('D7odd_cone.json'))['D_cone'])
L = core.parse('X0+2*X1+3*X2'); M = core.parse('Y0+3*Y1+5*Y2')
out = {}
for p in [2147483647, 2147483587]:
    for name, P in [('L*D7', mul(L, D7)), ('M*D7odd', mul(M, D7o)), ('X0^3*D7', mul(core.parse('X0**3'), D7)), ('Y0^3*D7odd', mul(core.parse('Y0**3'), D7o)),
                    ('X0^2*D7', mul(core.parse('X0**2'), D7)), ('A*D7', mul(A, D7)), ('A*D7odd', mul(A, D7o)), ('M*D7', mul(M, D7))]:
        ok, info = jacobian.in_jacobian(F, P, p)
        out[f'{name}@{p}'] = bool(ok)
        print(f'p={p}: {name} in J: {ok}', flush=True)
json.dump(out, open('verify2.json', 'w'))
