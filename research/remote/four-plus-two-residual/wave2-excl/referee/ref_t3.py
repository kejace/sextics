"""T3 referee checks: (a) character classes of the pairs (a,c) under G (order 24) -> number of eigenvector Jacobian families;
(b) Singular minAssGTZ / dim on the agent's saved annihilator quadrics (10-param and 8-param families)."""
import sys, json, subprocess, re, itertools; sys.path.insert(0, '/tmp/hodge/work/referee42')
import reflib as R
LOG = '/tmp/hodge/work/referee42/out/ref_t3.jsonl'
# generators as weight vectors mod 6 on (X0,X1,X2,Y0,Y1,Y2): h=(0,0,1,0,0,1), r2=(0,0,3,0,0,0), s=(0,0,0,3,3,0)
gens = [(0,0,1,0,0,1), (0,0,3,0,0,0), (0,0,0,3,3,0)]
def wt(e): return tuple(sum(a*b for a,b in zip(e,g)) % 6 for g in gens)
F, A = R.load_instance()
# check F is an eigenvector (all monomials same weight) and record variable weights
assert len({wt(e) for e in F}) == 1, 'F not G-eigenvector'
zw = [wt(tuple(1 if j == i else 0 for j in range(6))) for i in range(6)]
diffs = {}
for a in range(6):
    for c in range(6):
        d = tuple((zw[c][g] - zw[a][g]) % 6 for g in range(3))
        diffs.setdefault(d, []).append((a, c))
classes = sorted(diffs.items(), key=lambda t: -len(t[1]))
R.log(LOG, {'exp': 'eigenvector-classes', 'variable_weights(h,r2,s)': zw, 'n_classes': len(classes), 'class_sizes': [len(v) for _, v in classes],
            'classes': [(k, v) for k, v in classes], 'total_pairs': sum(len(v) for _, v in classes)})
# (b) Singular on saved quadrics
base = '/tmp/hodge/results/four-plus-two-residual/wave2-excl/receipts/'
for tag, fn, cmd in (('10-param chiF', 'excl-d6jac2-local-relax-chiF__classes__X0_X1___Y0_Y1___X2___Y2___10_params__contains.json', 'minass'),
                     ('8-param Y-odd merged', 'excl-d6jac2-local-relax-MERGED_Y_odd___Y0_Y1__x_grad__X0_X1_____X0_X1__x_grad__Y0_Y1.json', 'dim')):
    J = json.load(open(base + fn)); nP = len(J['pairs']); quads = J['quadrics']
    qs = [re.sub(r'\*\*', '^', q) for q in quads]
    vars_ = ','.join(f'w{i}' for i in range(nP))
    if cmd == 'minass':
        s = f'LIB "primdec.lib"; ring r=0,({vars_}),dp; ideal I={",".join(qs)}; ideal G=std(I); print(dim(G)); list L=minAssGTZ(I); print(size(L)); print(L); quit;'
    else:
        s = f'ring r=0,({vars_}),dp; ideal I={",".join(qs)}; ideal G=std(I); print(dim(G)); print(lead(G)); quit;'
    path = f'/tmp/hodge/work/referee42/out/ref_t3_{cmd}.sing'; open(path, 'w').write(s)
    out = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=1800).stdout
    R.log(LOG, {'exp': 'rank-one-locus-Singular', 'family': tag, 'n_params': nP, 'dim_W': J['dim_W'], 'n_quadrics': len(quads), 'singular': out.strip().replace('\n', ' | ')[:900]})
