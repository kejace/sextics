"""Slice a relaxation system by the spurious subspace F*S_{d-6}.

The relaxation quadrics q_r(a) = a^T S_r a (from relax.py, family order: J-type generators m F_i for i in 0..5,
m in monomials(d-5, (pD - par_i) % 2), then kernel reps) are invariant under a -> a + s for s in the coordinates
of F*G (G in S_{d-6}^{pD}) with lift A G Z/6.  This script (1) verifies S_r s = 0 mod p for every spurious vector,
(2) chooses one triangular coordinate per spurious vector and sets it to zero (a complement of the spurious
subspace), (3) writes the sliced msolve file.
  python3 wave2/slice_ms.py --d 6 --pD 0 --prime 1073741789 --ms out/relax_d6_p0_q.ms --out out/relax_d6_p0_q_slice.ms
"""
import sys, os, re, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *

ap = argparse.ArgumentParser()
ap.add_argument('--d', type=int, required=True)
ap.add_argument('--pD', type=int, default=0)
ap.add_argument('--prime', type=int, required=True)
ap.add_argument('--ms', required=True)
ap.add_argument('--out', required=True)
args = ap.parse_args()
F, A, q = load()
d, pD, p = args.d, args.pD, args.prime
Fp = core.reduce_mod(F, p)
Fi = [pderiv(Fp, i, p) for i in range(NV)]
# family order as in relax.py (J-type part); kernel reps come after and have zero spurious coordinates
fam = []
for i in range(NV):
    par = (pD - par_idx(i)) % 2
    for m in monomials(d - 5, par):
        fam.append((i, m))
nJ = len(fam)
index = {(i, m): k for k, (i, m) in enumerate(fam)}
lines = open(args.ms).read().split('\n')
names = lines[0].split(','); K = len(names)
assert int(lines[1]) == p
body = '\n'.join(lines[2:]).strip().rstrip(',')
polys = [s.strip() for s in body.split(',\n') if s.strip()]
# spurious vectors: F*G = (1/6) sum_i Z_i G F_i  -> coordinate 1/6 * coeff at generator (Z_i * mono of G) F_i
inv6 = pow(6, -1, p)
SG = monomials(d - 6, (pD - core.parity_of(F)) % 2) if d >= 6 else []
spur = []
for g in SG:
    v = [0] * K
    for i in range(NV):
        m = list(g); m[i] += 1; m = tuple(m)
        v[index[(i, m)]] = (v[index[(i, m)]] + inv6) % p
    spur.append(v)
print(f'[slice] d={d} pD={pD} K={K} (J-type {nJ}) spurious vectors {len(spur)}', file=sys.stderr)
# parse quadrics into coefficient dicts {(m,n): c}
term_re = re.compile(r'^(-?\d+)\*a(\d+)(?:\^2|\*a(\d+))?$')
quads = []
for s in polys:
    s = s.replace('\n', '').replace('-', '+-')
    coeffs = {}
    for t in s.split('+'):
        t = t.strip()
        if not t:
            continue
        mobj = term_re.match(t)
        assert mobj, t
        c = int(mobj.group(1)) % p; m = int(mobj.group(2)); n = int(mobj.group(3)) if mobj.group(3) else m
        if '^2' not in t and mobj.group(3) is None:
            raise ValueError(t)
        coeffs[(min(m, n), max(m, n))] = (coeffs.get((min(m, n), max(m, n)), 0) + c) % p
    quads.append(coeffs)
# invariance check: for each quadric, the bilinear form S_r(s, e_k) = 0 for all k, each spurious s
bad = 0
for r_, Qc in enumerate(quads):
    for v in spur:
        # (S v)_k = sum_n S[k,n] v_n with S symmetric: S[m,m] = c_mm, S[m,n] = c_mn/2
        Sv = [0] * K
        for (m, n), c in Qc.items():
            if m == n:
                Sv[m] = (Sv[m] + c * v[m]) % p
            else:
                h = c * pow(2, -1, p) % p
                Sv[m] = (Sv[m] + h * v[n]) % p
                Sv[n] = (Sv[n] + h * v[m]) % p
        if any(Sv):
            bad += 1
print(f'[slice] invariance violations: {bad} of {len(quads) * len(spur)} (quadric, spurious) pairs', file=sys.stderr)
# triangular slice coordinates: for spurious vector of G = monomial g, use generator (X0 * g) F_X0 (i = 0)
drop = []
for g in SG:
    m = list(g); m[0] += 1
    drop.append(index[(0, tuple(m))])
print(f'[slice] dropping coordinates {drop} ({[names[k] for k in drop]})', file=sys.stderr)
keep = [k for k in range(K) if k not in set(drop)]
newname = {k: f'b{t}' for t, k in enumerate(keep)}
out_polys = []
for Qc in quads:
    terms = []
    for (m, n), c in Qc.items():
        if m in drop or n in drop:
            continue
        terms.append(f'{c}*{newname[m]}^2' if m == n else f'{c}*{newname[m]}*{newname[n]}')
    if terms:
        out_polys.append('+'.join(terms))
with open(args.out, 'w') as fh:
    fh.write(','.join(newname[k] for k in keep) + '\n' + f'{p}\n' + ',\n'.join(out_polys) + '\n')
print(json.dumps({'batch': 'slice', 'd': d, 'pD': pD, 'p': p, 'K': K, 'n_spurious': len(spur), 'invariance_violations': bad,
                  'dropped': drop, 'n_quadrics': len(out_polys), 'vars_after': len(keep), 'out': args.out}))
