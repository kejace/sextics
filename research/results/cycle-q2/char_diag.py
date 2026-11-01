"""Character diagnostic for cone-type supports: which (2,0)-forms on exceptional divisors carry the
target character?  Condition: e_A - e_form - e_D has zero weight under the relevant group."""
import sys, itertools
sys.path.insert(0, '/tmp/hodge/results/cycle-q2')
import eq, numpy as np
from hodge import core
F, A, q, rec = core.load_instance('cycle-q2')
Gw = eq.symmetry_group(F); gens = eq.generators(Gw)
eA = np.array([2, 2, 3, 2, 2, 1])
def key(e):
    return int(eq.keys_of(np.array([e], dtype=np.int64), gens)[0])
def exps(deg, vars_):
    return [e for e in core.monomials(deg) if all(e[i] == 0 for i in range(6) if i not in vars_)]
X = [0, 1, 2]; Y = [3, 4, 5]
def vertex_plane_count(eD, deg, base_vars):
    other = Y if base_vars == X else X
    hits = []
    for n in exps(deg - 3, base_vars):
        for m in exps(3, other):
            form = np.array(n) + np.array(m)
            if np.all(((eA - form - np.array(eD)) @ Gw.T) % 6 == 0):
                hits.append((core.to_str({tuple(n): 1}), core.to_str({tuple(m): 1})))
    return hits, len(exps(deg - 3, base_vars)) * len(exps(3, other))
for name, eD, deg, bv in [('X-block sum X_i^6 (base: Fermat sextic curve in P^2_X, curve C_Y)', (6,0,0,0,0,0), 6, X),
                          ('Y-block sum Y_i^6', (0,0,0,6,0,0), 6, Y)]:
    hits, tot = vertex_plane_count(eD, deg, bv)
    print(f'{name}: forms on C x B = {tot}, with target character: {len(hits)} {hits[:6]}')
for deg in [3, 4, 5, 6]:
    for base_vars, nm in [(X, 'X'), (Y, 'Y')]:
        mons = exps(deg, base_vars)
        classes = {}
        for e in mons:
            classes.setdefault(key(e), []).append(e)
        good = []
        for k_, cl in classes.items():
            hits, tot = vertex_plane_count(cl[0], deg, base_vars)
            if hits:
                good.append((core.to_str({c: 1 for c in cl}), len(hits)))
        print(f'deg-{deg} cones in {nm} ({len(mons)} monomials, {len(classes)} classes): classes admitting target-character forms on C x B: {len(good)} -> {good[:8]}')
def vertex_line_hits(eD, deg, Dvars):
    line = [i for i in range(6) if i not in Dvars]; u, v = line
    HP = Gw[(Gw[:, u] - Gw[:, v]) % 6 == 0]
    eline = np.zeros(6, dtype=int); eline[u] = 1; eline[v] = 1
    hits = [core.to_str({tuple(n): 1}) for n in exps(deg - 4, Dvars)
            if np.all((((eA + eline) - np.array(n) - np.array(eD)) @ HP.T) % 6 == 0)]
    return hits, len(exps(deg - 4, Dvars)), len(HP)
for quad in itertools.combinations(range(6), 4):
    eD = tuple(6 if i == quad[0] else 0 for i in range(6))
    hits, tot, h = vertex_line_hits(eD, 6, list(quad))
    print(f'sextic Fermat cone on {",".join(core.VARS[j] for j in quad)} (|H_P|={h}): forms/point {tot}, target-character: {len(hits)} {hits[:6]}')
for deg in [4, 5, 6]:
    for quad in itertools.combinations(range(6), 4):
        mons = exps(deg, list(quad))
        classes = {}
        for e_ in mons:
            classes.setdefault(key(e_), []).append(e_)
        good = []
        for k_, cl in classes.items():
            hits, tot, h = vertex_line_hits(cl[0], deg, list(quad))
            if hits:
                good.append((core.to_str({c: 1 for c in cl})[:60], len(hits)))
        print(f'deg-{deg} cones on {",".join(core.VARS[j] for j in quad)}: {len(classes)} classes; with target-character forms: {len(good)} -> {good[:4]}')
