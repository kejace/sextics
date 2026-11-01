"""Batch 1 on cyclic-triple-residual: structural pieces, Fermat-type sums, cubic cones over all triples."""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import ctr
from ctr import *
out = []
def run(D, fam, note='', Dexpr=None, mode='auto'):
    r = test_D(D, fam, note=note, Dexpr=Dexpr, mode=mode)
    print(fmt(r), flush=True); out.append(r); return r
pc = pieces()
# --- structural pieces (degree <= 6) ---
run(pc['C'], 'struct:C', 'cubic coefficient C of Z in F = Z^2 + tCZ + H', Dexpr='C')
run(pc['U'], 'struct:U', 'U = X2^3 + (t/2)C, F = U^2 + Hp', Dexpr='U=X2**3+(5/14)*C')
run(pc['H'], 'struct:H', 'H = F|_{X2=0}', Dexpr='H')
run(pc['Hp'], 'struct:Hp', "H' = H - t^2C^2/4, F = U^2 + H'", Dexpr="H'")
run(pc['Delta'], 'struct:Delta', 'discriminant t^2C^2-4H', Dexpr='Delta')
run('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6', 'struct:fermat-part', 'pure-power part of F')
run('X0**6 + X1**6 + X2**6', 'struct:sumX')
run('Y0**6 + Y1**6 + Y2**6', 'struct:sumY')
run('X2**6 + Y0**6 + Y1**6 + Y2**6', 'struct:fqlx-reduced-support', 'fqlx excluded over C (1248 cols)')
run('X0**6 + Y0**6', 'struct:X0Y0-sum6', 'Fermat cone hyperplanes analogue')
# F minus each coupling term, and the coupling part alone
coup = ['5*X0**3*X1*Y0*Y1/7', '5*X0**3*X2**3/7', '5*X0**2*Y0*Y1*Y2**2/7', '5*X0*X1*Y1**3*Y2/7', '5*X0*X2**3*Y1*Y2/7', '5*X1*X2**3*Y0*Y1/7']
run(' + '.join(coup), 'struct:couplings', 'sum of the six coupling terms')
for c in coup:
    run(core.sub(F, P(c)), 'struct:F-minus-coupling', 'F minus ' + c, Dexpr='F - (' + c + ')')
# block radial derivatives (X_i d_i + Y_i d_i)F and partials
for i in range(3):
    Xi = [0]*6; Xi[i] = 1
    R = core.add(core.mul(P(V[i]), core.pderiv(F, i)), core.mul(P(V[i+3]), core.pderiv(F, i+3)))
    run(R, 'struct:block-radial', f'(X{i} d_X{i} + Y{i} d_Y{i})F', Dexpr=f'block-radial-{i}')
for i in range(6):
    run(core.pderiv(F, i), 'struct:partial', f'F_{V[i]} (degree 5, in J)', Dexpr=f'dF/d{V[i]}')
# --- Fermat-type sums of sixth powers over subsets, coefficient 1 ---
for k in (2, 3, 4):
    for S in itertools.combinations(range(6), k):
        Dexpr = ' + '.join(f'{V[i]}**6' for i in S)
        run(Dexpr, f'sum6-{k}', f'subset {[V[i] for i in S]}')
# --- cubic cones over all triples (definite parity) ---
def cubic_cones(S):
    names = [V[i] for i in S]
    Xs = [n for n in names if n[0] == 'X']; Ys = [n for n in names if n[0] == 'Y']
    cands = []
    if len(Ys) == 0 or len(Ys) == 3:
        cands.append(' + '.join(f'{n}**3' for n in names))
        cands.append(f'{names[0]}**3 + {names[1]}**3 + {names[2]}**3 + {names[0]}*{names[1]}*{names[2]}')
    if len(Ys) == 1:
        a, b = Xs; y = Ys[0]
        cands.append(f'{a}**3 + {b}**3 + {a}*{y}**2')
        cands.append(f'{a}**3 + {b}**3 + {b}*{y}**2')
        cands.append(f'{a}**3 + {b}**3 + ({a}+{b})*{y}**2')
    if len(Ys) == 2:
        x = Xs[0]; a, b = Ys
        cands.append(f'{a}**3 + {b}**3 + {x}**2*{a}')
        cands.append(f'{a}**3 + {b}**3 + {x}**2*{b}')
        cands.append(f'{a}**3 + {b}**3 + {x}**2*({a}+{b})')
    return cands
for S in itertools.combinations(range(6), 3):
    for Dexpr in cubic_cones(S):
        run(Dexpr, 'cone3', f'cubic cone over {[V[i] for i in S]}')
json.dump(out, open('batch1.json', 'w'), indent=1)
print('BATCH1 DONE')
