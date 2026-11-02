"""Independent check of an msolve system built by fpt_msolve.py for a control instance: substitute the (gauge-fixed)
known witness a_w mod p and verify the equations become a consistent linear system in x (mod p)."""
import sys, re
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-search')
from fptlib import *
from fpt_msolve import modp_family, nm, table, particular_solutions
from flint import nmod_mat
iid, msfile, p = sys.argv[1], sys.argv[2], int(sys.argv[3])
inst0 = Instance(iid); d, Dw, Nw, Qw = core.load_witness(iid)
subg = symmetry.stabilizer(inst0.G, [inst0.A, Dw]); inst = Instance(iid, gens=subg)
lab = label_of_poly(inst, Dw)
Dm, Bt = modp_family(inst, d, lab, p); K0, K = len(Dm), len(Bt)
didx = {m: j for j, m in enumerate(Dm)}
Dp = reduce_mod(Dw, p)
cvec = [[Dp.get(m, 0) % p] for m in Dm]
B = [[Bt[k][j] for k in range(K)] for j in range(K0)]
aw = particular_solutions(B, cvec, p)
assert aw is not None, 'witness not in the mod-p family'
aw = [row[0] for row in aw]
lines = open(msfile).read().strip().split('\n')
varnames = lines[0].split(','); assert int(lines[1]) == p
eqs = [l.rstrip(',') for l in lines[2:] if l.strip()]
nC = sum(1 for v in varnames if v.startswith('x'))
# F-gauge: find lambda with s.(a_w + lam * s_dir) = 0 is implicit: instead solve the gauge/normalisation equations for the
# representative: we test consistency of ALL equations treating a = t*a_w + u*s (s = F-multiple direction) is complex; so
# we simply substitute a = a_w and drop the F-gauge/normalisation linear equations in a, checking the remaining ones.
def parse_term(t):
    parts = t.split('*'); coef = int(parts[0]); vs = []
    for q in parts[1:]:
        if '^' in q:
            v, e = q.split('^'); vs += [v] * int(e)
        else:
            vs.append(q)
    return coef, vs
rows = []; rhs = []; n_dropped = 0; n_quad = 0
for eq in eqs:
    terms = [parse_term(t) for t in eq.split('+')]
    only_a_linear = all(len(vs) <= 1 and all(v.startswith('a') for v in vs) for c, vs in terms)
    if only_a_linear:
        n_dropped += 1; continue
    row = [0] * nC; const = 0
    for c, vs in terms:
        val = c
        xs = [v for v in vs if v.startswith('x')]
        for v in vs:
            if v.startswith('a'):
                val = val * aw[int(v[1:])] % p
        if len(xs) == 0:
            const = (const + val) % p
        elif len(xs) == 1:
            row[int(xs[0][1:])] = (row[int(xs[0][1:])] + val) % p
        else:
            raise ValueError('quadratic in x?')
    rows.append(row); rhs.append([(-const) % p])
sol = particular_solutions(rows, rhs, p) if rows else []
print(f'{iid}: K={K} nC={nC} equations={len(eqs)} (dropped {n_dropped} linear-in-a gauge/normalisation eqs); '
      f'linear system in x at a_w: {len(rows)} eqs -> {"CONSISTENT (witness satisfies the msolve system)" if sol is not None else "INCONSISTENT (bug!)"}')
