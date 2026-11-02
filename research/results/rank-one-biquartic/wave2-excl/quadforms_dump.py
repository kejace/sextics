"""Dump (mod p) the quadratic forms in a that lie in the J-relaxed weight-2 ideal (both parity blocks) for degree d,
as a Singular script that tests whether their common zero locus is contained in the spurious subspace Sp = F*S_{d-6}.
Usage: python3 quadforms_dump.py <d> <prime> [--nonj FILE]"""
import sys, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/rank-one-biquartic/wave2-excl')
from w2lib import *
from w2log import log
ap = argparse.ArgumentParser(); ap.add_argument('d', type=int); ap.add_argument('prime', type=int); ap.add_argument('--nonj', default=''); ap.add_argument('--tag', default='')
a = ap.parse_args(); d, p = a.d, a.prime
tag = a.tag or f'quadforms_d{d}_p{p}'
basis = polar_basis(d)
if a.nonj:
    for k, rec_ in enumerate(json.load(open(a.nonj))):
        basis.append((('K', rec_.get('label', k)), core.parse(rec_['D']), [core.parse(s) for s in rec_['N']]))
K = len(basis); labidx = {b[0]: k for k, b in enumerate(basis)}
HF = HessianFamily(d, basis)
forms = {0: [], 1: []}
ranks = {}
for blk in (0, 1):
    info, Mr, piv = HF.analyze(blk, p)
    B = HF.blocks[blk]; C = B['ncols']; nJ = HF.nJ[blk]
    for i, j in enumerate(piv):
        if j < nJ: continue
        row = {k: int(Mr[i, k]) for k in range(j, C) if int(Mr[i, k])}
        assert all(k >= nJ for k in row)
        forms[blk].append({B['cols'][k][1:]: v for k, v in row.items()})
    ranks[blk] = dict(rank=info['rank'], rank_J=info['rank_J'], n_forms=len(forms[blk]))
    del Mr
print(json.dumps(ranks), flush=True)
# spurious groups
sp_groups = []
if d >= 6:
    for m in monomials(d - 6):
        grp = [labidx[('J', i, tuple(m[j] + (1 if j == i else 0) for j in range(NV)))] for i in range(NV)]
        sp_groups.append(grp)
in_group = {k for grp in sp_groups for k in grp}
coords = [(f'a{k}', {k: 1}) for k in range(K) if k not in in_group]
for grp in sp_groups:
    for i in range(1, NV):
        coords.append((f'a{grp[0]}-a{grp[i]}', {grp[0]: 1, grp[i]: -1}))
# Singular script: ideal I of all forms; for each coordinate u: dim(std(I + (u - 1))) == -1 ?  (weighted-homogeneous cone => u in radical)
lines = [f'ring r = {p}, (' + ','.join(f'a{k}' for k in range(K)) + '), dp;', 'ideal I;']
for blk in (0, 1):
    for fm in forms[blk]:
        terms = []
        for (k, l), v in fm.items():
            terms.append(f'{v}*a{k}*a{l}' if k != l else f'{v}*a{k}^2')
        lines.append('I[size(I)+1] = ' + '+'.join(terms) + ';')
lines.append('ideal G = std(I);')
lines.append('"dim V(I) (affine):", dim(G), " #forms:", size(I);')
lines.append('int bad = 0; int k; ideal T;')
for name, lin in coords:
    u = '+'.join(f'({v})*a{k}' for k, v in lin.items())
    lines.append(f'T = G, {u} - 1; if (dim(std(T)) >= 0) {{ bad++; "coordinate {name} NOT in radical"; }}')
lines.append('"coordinates not in radical:", bad, "of", ' + str(len(coords)) + ';')
lines.append('quit;')
open(f'{PK}/{tag}.sing', 'w').write('\n'.join(lines) + '\n')
json.dump(dict(ranks=ranks, K=K, coords=[n for n, _ in coords], forms={str(b): [[list(k), v] for k, v in fm.items()] for b in (0, 1) for fm in [dict((tuple(kk), vv) for kk, vv in f_.items()) for f_ in forms[b]]}), open(f'{PK}/{tag}.json', 'w'))
log(dict(exp='hessian-relaxation-quadforms', d=d, prime=p, K=K, forms_blk0=ranks[0]['n_forms'], forms_blk1=ranks[1]['n_forms'], rank_blk0=ranks[0]['rank'], rank_blk1=ranks[1]['rank'], job=tag,
         note='quadratic forms in a (mod p) spanning the J-relaxed weight-2 ideal; Singular script tests V(forms) subset Sp'), path=f'{PK}/{tag}.log.jsonl')
print('wrote', f'{PK}/{tag}.sing')
