"""Ages (Reid-Tai) of every element of G = <sigma, tau> (Z/14) at every fixed component on X.
Element g = (sigma^k, tau^j) acts by weights (k, 2k, 4k, 0,0,0)/7 + j*(1/2,1/2,1/2,0,0,0) (tau realised as X -> -X).
Fixed components: coordinate points e_x, e_y, e_z (tangent space = kernel of dF), the curve {X=0, G(Y)=0}
(normal directions x,y,z), and the curve {Y=0, K=0} (normal directions t,r,s; fixed by tau only).
Age at a point = sum of fractional parts of the weights on the tangent space of X in the affine chart (weights relative
to the chart coordinate)."""
from fractions import Fraction as Fr
w7 = [1,2,4,0,0,0]
def weights(k, j):
    return [Fr(k*w7[i], 7) + (Fr(j,2) if i < 3 else 0) for i in range(6)]
def frac(a): return a - (a.numerator // a.denominator)
out = []
# coordinate points: chart index c, tangent directions = all other indices except the one killed by dF
tangent = {0: [2,3,4,5], 1: [0,3,4,5], 2: [1,3,4,5]}   # at e_x: dF = x^5 dy -> drop y; e_y: drop z; e_z: drop x
for k in range(7):
    for j in range(2):
        if k == 0 and j == 0: continue
        wt = weights(k, j)
        row = {'g': f'sigma^{k} tau^{j}'}
        for c, tdirs in tangent.items():
            age = sum(frac(wt[i] - wt[c]) for i in tdirs)
            row[f'age@e_{"xyz"[c]}'] = str(age)
        # curve {X=0, G=0}: fixed by all g; normal weights relative to a Y-chart coordinate (weight 0): x,y,z
        row['age@{X=0,G=0}'] = str(sum(frac(wt[i]) for i in range(3)))
        # curve {Y=0, K=0}: fixed only if g acts trivially on P^2_X projectively i.e. k == 0 (tau): normal weights t,r,s rel. x
        if k == 0:
            row['age@{Y=0,K=0}'] = str(sum(frac(wt[i] - wt[0]) for i in range(3,6)))
        out.append(row)
for r in out: print(r)
junior = [r['g'] for r in out if any(Fr(v) < 1 for key, v in r.items() if key.startswith('age'))]
print('JUNIOR elements (age < 1 somewhere):', junior)
