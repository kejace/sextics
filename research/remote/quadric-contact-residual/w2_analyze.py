#!/usr/bin/env python3
"""Inspect the best minima saved by w2_varpro.py / w2_hess.py (varpro-*.json, hess-*.json) in the remote mirror:
which family members dominate, how the best D looks (largest monomial coefficients), and whether minima repeat.
  python3 w2_analyze.py [pattern]
"""
import sys, os, json, glob
import numpy as np
MIR = '/tmp/hodge/remote/quadric-contact-residual'
pat = sys.argv[1] if len(sys.argv) > 1 else 'varpro-*.json'
for path in sorted(glob.glob(os.path.join(MIR, pat))):
    J = json.load(open(path))
    rec, members, best = J['rec'], J['members'], J['best']
    print(f'== {os.path.basename(path)}: d={rec["d"]} chi={rec["chi"]} K={rec.get("K", rec.get("n_members"))} dimP={rec.get("dimP")} '
          f'restarts={rec.get("restarts_done")} best={rec.get("best_rel", rec.get("best_frel"))}')
    if not best:
        continue
    A = np.array([b['a'] for b in best])
    vals = [b.get('rel', b.get('frel')) for b in best]
    # cluster: pairwise |cos|
    G = np.abs(A @ A.T)
    print('   top values:', ', '.join(f'{v:.2e}' for v in vals[:8]))
    print('   |cos| between best and others:', ', '.join(f'{G[0, j]:.2f}' for j in range(1, min(8, len(best)))))
    a = A[0]
    idx = np.argsort(-np.abs(a))[:8]
    print('   dominant members of the best minimum:', ', '.join(f'{members[i]}:{a[i]:+.3f}' for i in idx))
