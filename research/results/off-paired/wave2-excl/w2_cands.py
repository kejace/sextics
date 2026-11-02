"""Candidate-denominator lists for the wave-2 re-verification batches (expression strings over block names).
Writes JSONL files {"D": expr, "family": tag, "label": name} into the packet dir."""
import sys, json, itertools
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib
from hodge import core
from hodge.core import monomials, NV, VARS

PK = w2lib.PK
B = dict(w2lib.blocks()['_polys'])
for i, v in enumerate(VARS):
    B[v] = core.parse(v)
bchar = {k: w2lib.char_of(v) for k, v in B.items()}
bdeg = {k: core.degree(v) for k, v in B.items()}

def mono_str(m):
    return '*'.join(f'{VARS[i]}**{k}' if k > 1 else VARS[i] for i, k in enumerate(m) if k) or '1'

def write(path, cands):
    with open(path, 'w') as f:
        for c in cands:
            f.write(json.dumps(c) + '\n')
    print(path, len(cands))

# ---- wave-1 structured families (attack1 §2.3, 2.5, 2.6) -------------------------------------------
w1 = []
parts = ['F' + v for v in VARS]
for P in parts:
    for m in monomials(2):
        w1.append({'D': f'{P}*{mono_str(m)}', 'family': 'partial*monomial', 'label': f'{P}*{mono_str(m)}'})
for m in monomials(3):
    w1.append({'D': f'FX0*{mono_str(m)}', 'family': 'partial*monomial3', 'label': f'FX0*{mono_str(m)}'})
    w1.append({'D': f'D0*{mono_str(m)}', 'family': 'quartic-colon*monomial', 'label': f'D0*{mono_str(m)}'})
# subset Euler fields E_I(F) = sum_{i in I} Z_i F_i, |I| = 1..3, times a variable
for r in (1, 2, 3):
    for I in itertools.combinations(range(NV), r):
        E = '+'.join(f'{VARS[i]}*F{VARS[i]}' for i in I)
        for j in range(NV):
            w1.append({'D': f'({E})*{VARS[j]}', 'family': 'subset-euler*monomial', 'label': f'E{"".join(str(i) for i in I)}*{VARS[j]}'})
# Fermat-deformation blocks
for blk in ('f6', 'g6', '(f6+g6)', 'Qx'):
    for dm in (1, 2):
        for m in monomials(dm):
            w1.append({'D': f'{blk}*{mono_str(m)}', 'family': 'fermat-deformation*monomial', 'label': f'{blk}*{mono_str(m)}'})
for a, b in [('f6','g6'),('f6','Qx'),('g6','Qx'),('f6','H'),('g6','H'),('Qx','H'),('(f6+g6)','Qx'),('(f6+g6)','H'),('f6','(f6+g6)'),('g6','(f6+g6)'),
             ('f6','C'),('g6','C'),('(f6+g6)','C'),('Qx','C'),('f6','D0'),('g6','D0'),('(f6+g6)','D0'),('Qx','D0')]:
    w1.append({'D': f'{a}*{b}', 'family': 'fermat-deformation*block', 'label': f'{a}*{b}'})
w1.append({'D': '27*H**2-2*C**3', 'family': 'genus2-discriminant-component', 'label': '27H^2-2C^3'})
w1.append({'D': '27*H**2+2*C**3', 'family': 'genus2-discriminant-component', 'label': '27H^2+2C^3'})
w1.append({'D': 'f6', 'family': 'fermat-cone', 'label': 'f6'})
w1.append({'D': 'g6', 'family': 'fermat-cone', 'label': 'g6'})
w1.append({'D': 'f6+g6', 'family': 'fermat-cone', 'label': 'f6+g6'})
w1.append({'D': 'Qx', 'family': 'mixed-part', 'label': 'Qx'})
w1.append({'D': 'H', 'family': 'base-block', 'label': 'H'})
w1.append({'D': 'C', 'family': 'base-block', 'label': 'C'})
w1.append({'D': 'D0', 'family': 'colon-quartic', 'label': 'D0'})
w1.append({'D': 'X2**6-Y2**6', 'family': 'mirror-product', 'label': 'M'})
w1.append({'D': 'X2**6+Y2**6', 'family': 'mirror-product', 'label': 'S'})
# 2x2 Hessian minors (E1 only in wave 1)
Fp = B['F']
for (i, k), (j, l) in itertools.product(itertools.combinations(range(NV), 2), repeat=2):
    pass
hess = []
seen = set()
for i, j in itertools.combinations(range(NV), 2):
    for k, l in itertools.combinations(range(NV), 2):
        key = (i, j, k, l)
        if key in seen: continue
        seen.add(key)
        # minor rows i,j cols k,l: F_ik F_jl - F_il F_jk ; write via second derivatives of block expressions is awkward -> expand numerically
        from hodge.core import pderiv, mul, sub, to_str
        Fik = pderiv(pderiv(Fp, i), k); Fjl = pderiv(pderiv(Fp, j), l); Fil = pderiv(pderiv(Fp, i), l); Fjk = pderiv(pderiv(Fp, j), k)
        Mn = sub(mul(Fik, Fjl), mul(Fil, Fjk))
        if Mn:
            hess.append({'D': to_str(Mn), 'family': 'hessian-2x2-minor', 'label': f'hess[{i}{j}|{k}{l}]'})
write(f'{PK}/cands_wave1.jsonl', w1)
write(f'{PK}/cands_hessian.jsonl', hess)

# ---- dossier E2 block products (squarefree products of <= 3 blocks, semi-invariant, degrees 7-14) ----
blk = ['X0','X1','X2','Y0','Y1','Y2','P2','M','S','M2','S2','D0','H','C','R'] + parts
cands = []
seen = set()
for r in (1, 2, 3):
    for combo in itertools.combinations(blk, r):
        if sum(1 for b in combo if b.startswith('F')) > 1:
            continue
        if any(bchar[b] is None for b in combo):
            continue
        d = sum(bdeg[b] for b in combo)
        if not 7 <= d <= 14:
            continue
        expr = '*'.join(combo)
        if expr in seen: continue
        seen.add(expr)
        cands.append({'D': expr, 'family': 'block-product', 'label': expr, 'd': d})
from collections import Counter
print('block products by degree:', sorted(Counter(c['d'] for c in cands).items()))
write(f'{PK}/cands_blocks.jsonl', cands)
