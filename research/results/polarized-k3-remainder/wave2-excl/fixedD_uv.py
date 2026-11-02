"""Fixed-D driver in the UV frame: legality, E1, sector-restricted E2 (mu_6 x parity), per-character components of
Ahat for eigenvector D, and an exact Farkas certificate when inconsistent.  D is given as an expression in
X0..Y2 standing for U0..V2 (UV frame).  Usage: python3 fixedD_uv.py '<label>' '<family>' '<expr>' [--no-farkas] [--components]"""
import sys, time, json, argparse
sys.path.insert(0, '.')
from k3lib import *
from farkas import build_exact, farkas, save_certificate
from logutil import log
from hodge import symmetry, jacobian, fixed_d

ap = argparse.ArgumentParser()
ap.add_argument('label'); ap.add_argument('family'); ap.add_argument('expr')
ap.add_argument('--no-farkas', action='store_true'); ap.add_argument('--components', action='store_true')
ap.add_argument('--prime', type=int, default=1000000007); ap.add_argument('--frame', default='uv')
args = ap.parse_args()
t0 = time.time()
if args.frame == 'uv':
    F, A, q = load_uv()
else:
    F, A, q = load_xy()
D = core.parse(args.expr)
d = core.degree(D)
assert core.is_homogeneous(D)
p = args.prime
# legality: F | D ?
divisible, _ = jacobian.divisible_by(F, D, p)
rec = dict(exp='fixedD', label=args.label, family=args.family, frame=args.frame, D=args.expr if len(args.expr) < 400 else args.expr[:400] + '...', d=d, terms=len(D), prime=p)
if divisible:
    rec.update(status='ILLEGAL (F | D mod p)'); log(**rec); sys.exit(0)
G = symmetry.diagonal_group(F)
chiD = [symmetry.character(D, w, G.order) for w in G.elements]
rec['h_eigenvector'] = all(c is not None for c in chiD)
rec['swap_parity'] = 'sym' if sub(swap(D), D) == {} else ('anti' if add(swap(D), D) == {} else 'none')
rec['y_parity'] = core.parity_of(D)
# E1
ok1, info1 = fixed_d.e1_solvable(F, A, q, D, p)
rec['e1'] = ok1; rec['e1_info'] = info1
if not ok1:
    rec.update(status='E1 fails (A D not in J(F) mod p)', consistent=False, seconds=round(time.time() - t0, 1)); log(**rec); print(rec); sys.exit(0)
# components of A to test: full A, and (if eigenvector and requested) the characters
targets = [('full', A)]
if args.components and rec['h_eigenvector']:
    for r, P in sorted(split_vdeg_mod(A).items()):
        targets.append((f'vdeg{r}', P))
results = {}
for name, At in targets:
    t1 = time.time()
    allowed, subg = symmetry.sector_filter(F, At, q, D, G)
    S = build_exact(F, At, q, D, allowed=allowed if subg else None, verbose=True)
    out = farkas(S, p=p, verbose=True) if not args.no_farkas else None
    if out is None:
        from farkas import dense_modp
        Mb = dense_modp(S, p); ra = Mb.rank(); r0 = dense_modp(S, p, with_rhs=False).rank()
        out = dict(prime=p, rows=S['nrows'], cols=S['ncols'], rank_modp=r0, rank_aug_modp=ra, status='consistent-modp' if ra == r0 else 'inconsistent-modp')
    res = {k: v for k, v in out.items() if k != 'certificate'}
    res.update(subgroup=len(subg), rows=S['nrows'], cols=S['ncols'], seconds=round(time.time() - t1, 1))
    if 'certificate' in out:
        path = f"cert-{args.label.replace(' ', '_').replace('/', '_')}-{name}.json"
        save_certificate(path, S, out, dict(label=args.label, family=args.family, frame=args.frame, D=args.expr, target=name))
        res['certificate_file'] = path
    results[name] = res
    print(name, res, flush=True)
rec.update(results=results, consistent=any(r['status'].startswith('consistent') for r in results.values()) if 'full' not in results else results['full']['status'].startswith('consistent'),
           status=results['full']['status'], seconds=round(time.time() - t0, 1))
log(**rec)
