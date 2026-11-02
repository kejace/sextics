"""Independent recomputation of the dossier's A.5 fixed-D exclusions (degree 5, 6, 10) with own code,
two new primes (999983, 1299709) and EXACT rank over Z (fmpz_mat) -> rigorous verdict."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/verify')
import vcheck as V
from fractions import Fraction
F, A, q = V.load_instance()
F2 = V.scale(F, 2)
X = lambda i: V.mono(**{f'X{i}': 1}); Y = lambda i: V.mono(**{f'Y{i}': 1})
Fi = [V.pderiv(F2, i) for i in range(6)]         # 2 F_i (integral)
cases = []
for i in range(3): cases.append((f'F_X{i} (x2)', Fi[i]))
for i in range(3): cases.append((f'F_Y{i} (x2)', Fi[3 + i]))
B6 = {}
for i in range(3): B6 = V.add(B6, V.mul(X(i), Fi[i]))
cases.append(('B6 = sum X_i F_X_i (x2)', B6))
for i in range(3):
    cases.append((f'G_rad,{i} = (X_i d_Xi + Y_i d_Yi)(2F)', V.add(V.mul(X(i), Fi[i]), V.mul(Y(i), Fi[3 + i]))))
for i in range(3):
    for j in range(3):
        cases.append((f'X{j} F_X{i} (x2)', V.mul(X(j), Fi[i])))
        cases.append((f'Y{j} F_Y{i} (x2)', V.mul(Y(j), Fi[3 + i])))
cases.append(('X0^2X1^2X2^2', V.mono(X0=2, X1=2, X2=2)))
cases.append(('Y0^2Y1^2Y2^2', V.mono(Y0=2, Y1=2, Y2=2)))
# degree 10
K10 = V.add(V.mul(V.mono(Y2=5), Fi[4]), V.mul(V.mono(Y1=5), Fi[5]), -1)
cases10 = [('K10 = (Y2^5 d_Y1 - Y1^5 d_Y2)(2F)', K10), ('F_X2^2 (x4)', V.mul(Fi[2], Fi[2])), ('F_X0 F_X1 (x4)', V.mul(Fi[0], Fi[1])),
           ('F_Y2^2 (x4)', V.mul(Fi[5], Fi[5])), ('F_Y0 F_Y1 (x4)', V.mul(Fi[3], Fi[4])), ('F_X2 F_X0 (x4)', V.mul(Fi[2], Fi[0]))]
which = sys.argv[1] if len(sys.argv) > 1 else 'low'
todo = cases if which == 'low' else cases10
for name, D in todo:
    D = {m: int(v) for m, v in D.items()}
    t0 = time.time()
    S = V.FixedD(F, A, q, D, verbose=True)
    r1 = S.ranks_mod(999983); r2 = S.ranks_mod(1299709)
    t1 = time.time()
    try:
        ex = S.ranks_exact(); texact = round(time.time() - t1, 1)
    except Exception as e:
        ex = ('ERR', str(e)); texact = None
    res = {'D': name, 'd': S.d, 'order': S.order, 'rows': S.nrows, 'cols': S.ncols, 'nN': S.nN, 'nQ': S.nQ,
           'rank_aug_999983': r1, 'rank_aug_1299709': r2, 'rank_aug_exact_Z': ex, 't_exact': texact, 'seconds': round(time.time() - t0, 1)}
    verdict = 'CONFIRMED-RIGOROUS' if (isinstance(ex[0], int) and ex[1] == ex[0] + 1) else ('CONFIRMED-MODP-ONLY' if r1[1] == r1[0] + 1 and r2[1] == r2[0] + 1 else 'REFUTED-OR-CONSISTENT')
    res['verdict'] = verdict
    print(json.dumps(res), flush=True)
    V.log(claim=f'fixed-D exclusion of D = {name} at q=2 (no certificate)', source='dossier A.5 (p=10007, 65521) / attack1 exp_g', verdict=verdict, evidence=res,
          notes='exact rank over Z via flint fmpz_mat: rank(aug) = rank(M)+1 proves inconsistency over Q hence C (R2); mod-p ranks at two new primes')
