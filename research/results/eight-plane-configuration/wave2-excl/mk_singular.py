import sys
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration/wave2-excl')
from w2lib import *
def sing(P):
    terms = []
    for e, c in sorted(P.items(), reverse=True):
        mon = '*'.join(f'{core.VARS[i]}^{k}' if k > 1 else core.VARS[i] for i, k in enumerate(e) if k)
        terms.append(f'({c})*{mon}' if mon else f'({c})')
    return ' + '.join(terms)
p = int(sys.argv[1])
print(f'''ring r = {p},(X0,X1,X2,Y0,Y1,Y2),dp;
option(redSB);
poly F = {sing(Fint)};
poly A = {sing(Aint)};
ideal J = jacob(F);
int t0 = timer;
ideal I = J, A;
ideal G = std(I);
print("std done, size " + string(size(G)) + " time " + string(timer - t0));
int k;
for (k = 0; k <= 17; k++) {{ print("HF(S/(J,A)) deg " + string(k) + " = " + string(size(kbase(G, k)))); }}
// Jacobian ring itself (check regular sequence mod p): Hilbert function of S/J
ideal GJ = std(J);
print("std(J) done, size " + string(size(GJ)) + " time " + string(timer - t0));
for (k = 0; k <= 25; k++) {{ print("HF(S/J) deg " + string(k) + " = " + string(size(kbase(GJ, k)))); }}
// quotient (J : A): its Hilbert function minus that of J gives dim Ann_R(A)_k
ideal Qt = quotient(GJ, ideal(A));
ideal GQ = std(Qt);
for (k = 0; k <= 12; k++) {{ print("dim Ann_R(A) deg " + string(k) + " = " + string(size(kbase(GJ, k)) - size(kbase(GQ, k)))); }}
quit;
''')
