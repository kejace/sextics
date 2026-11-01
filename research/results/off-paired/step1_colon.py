import oplib, time, json, sys
from hodge import core
from hodge.core import mul, monomials
F, A, q = oplib.load()
p = oplib.P0
chars = [(a, b) for a in range(6) for b in range(6)]
out = {}
for d in range(0, 11):
    t0 = time.time()
    row = {}
    tot = 0; totg = 0
    for c in chars:
        dim, basis, info = oplib.colon_space(F, A, q, d, c, p)
        # Euler gauge: subtract dim of F*S_{d-6} in character c (F has char (0,0)) -- these lie in colon space automatically
        g = len(oplib.monos_char(d - 6, c)) if d >= 6 else 0
        row[c] = (dim, dim - g, info['nD'])
        tot += dim; totg += dim - g
    out[d] = row
    nz = {f'{c[0]}{c[1]}': v for c, v in row.items() if v[0]}
    print(f'd={d}: total colon dim {tot}, after Euler gauge {totg}, time {time.time()-t0:.1f}s')
    print('   per char (dim, dim-gauge, #monomials):', nz)
    sys.stdout.flush()
json.dump({str(d): {f'{c[0]},{c[1]}': v for c, v in row.items()} for d, row in out.items()}, open('colon_dims.json', 'w'), indent=1)
