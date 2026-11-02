"""Idempotent: engine v3 -- rank/nullity of the COEFFICIENT matrix for inconsistent systems (dense + compressed)."""
import sys
P = '/home/kejace/devel/sextics-work/hodge/fixed_d.py'
f = open(P).read()
if "'rank_aug': rank_aug" in f and 'ENGINE_VERSION = 3' in f:
    print('fixed_d rankfix: already patched'); sys.exit(0)
if 'import os as _os' not in f:
    f = f.replace('from __future__ import annotations\nimport sys, time\n', 'from __future__ import annotations\nimport os as _os\nimport sys, time\n', 1)
old = '''        consistent = not (pivots and pivots[-1] == C)
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        res = {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'x': x,
               'pivots': pivots, 'rows': self.nrows, 'cols': C,
               't_build': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}'''
new = '''        consistent = not (pivots and pivots[-1] == C)
        rank_aug = rank
        if not consistent:
            rank -= 1            # rank of the coefficient matrix M (the last pivot is the rhs column)
            pivots = pivots[:-1]
        x = None
        if consistent:
            x = [0] * C
            for i, j in enumerate(pivots):
                x[j] = int(Mr[i, C])
        res = {'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'nullity': C - rank, 'x': x,
               'pivots': pivots, 'rows': self.nrows, 'cols': C,
               't_build': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}'''
if old in f: f = f.replace(old, new)
else: assert "'rank_aug': rank_aug, 'nullity': C - rank, 'x': x,\n               'pivots'" in f, 'dense solve block not found'
old = '''    consistent = not (pivots and pivots[-1] == C)
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots):
            x[j] = int(Mr[i, C])
    if verbose:
        print(f'[compressed p={S.p}] k={k} consistent={consistent} rank={rank} nullity={C-rank} '
              f'(compress {t1-t0:.1f}s, rref {t2-t1:.1f}s)', file=sys.stderr)
    return {'consistent': consistent, 'rank': rank, 'nullity': C - rank, 'x': x, 'pivots': pivots,
            'rows': S.nrows, 'cols': C, 'k': k, 't_compress': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}'''
new = '''    consistent = not (pivots and pivots[-1] == C)
    rank_aug = rank
    if not consistent:
        rank -= 1                # rank of the coefficient matrix M (the last pivot is the rhs column)
        pivots = pivots[:-1]
    x = None
    if consistent:
        x = [0] * C
        for i, j in enumerate(pivots):
            x[j] = int(Mr[i, C])
    if verbose:
        print(f'[compressed p={S.p}] k={k} consistent={consistent} rank={rank} nullity={C-rank} '
              f'(compress {t1-t0:.1f}s, rref {t2-t1:.1f}s)', file=sys.stderr)
    return {'consistent': consistent, 'rank': rank, 'rank_aug': rank_aug, 'nullity': C - rank, 'x': x, 'pivots': pivots,
            'rows': S.nrows, 'cols': C, 'k': k, 't_compress': round(t1 - t0, 2), 't_rref': round(t2 - t1, 2)}'''
if old in f: f = f.replace(old, new)
else: assert "'rank_aug': rank_aug, 'nullity': C - rank, 'x': x, 'pivots': pivots,\n            'rows'" in f, 'compressed block not found'
old = '''ENGINE_VERSION = 2
import sys as _sys
_sys.stderr.write('[hodge engine v2] NOTE: a sign bug in the F*N(D) term of FixedDSystem/exact solve was fixed on 2026-11-01 ~17:25 local time. '
                  'Any consistent/inconsistent/rank result computed with the engine BEFORE that time is unreliable and must be recomputed. '
                  'Record "engine": 2 in your logs from now on.\\n')'''
new = '''ENGINE_VERSION = 3
import sys as _sys
if not _os.environ.get('HODGE_QUIET'):
    _sys.stderr.write('[hodge engine v3] NOTE: engine v2 (2026-11-01 17:25) fixed a sign bug in the F*N(D) term; results from before '
                      'that time must be recomputed.  Engine v3 (2026-11-02 00:55 UTC) reports rank/nullity of the COEFFICIENT matrix '
                      'for inconsistent systems (v2 reported the augmented rank, i.e. rank+1 / nullity-1; consistency verdicts and all '
                      'consistent-case numbers were already correct).  Record "engine": 3 in your logs from now on.\\n')'''
if old in f: f = f.replace(old, new)
else: assert 'ENGINE_VERSION = 3' in f, 'version block not found'
open(P, 'w').write(f); print('fixed_d rankfix: patched')
