"""d=7 structural families (sigma- and 3-cycle-adapted), tested with the compressed target-image test."""
import sys, json, random, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
from families import family
from symproj import invariant_part, sign_part, t01_invariant, t01_anti, span_basis
from hodge.core import add, scale, parse, to_str
from fractions import Fraction as Fr
H = Fr(1,2)
Msig = [[H,H,0],[H,-H,0],[0,0,1]]
Mcyc = [[1,1,0],[1,0,1],[1,-1,-1]]
random.seed(11)
def rnd(polys, k=3):
    D = {}
    for P in polys: D = add(D, scale(P, Fr(random.randint(-k, k))))
    return D
def run(name, D, fam):
    e = test_D_image(name, fam, D=D)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'served', e.get('served_dim'), e['seconds'], 's', flush=True)
specs = [
 ('cyc (r1,r2,s1,s2)^4 d7 [quadruple pts at S3-fixed pts]', Mcyc, [1,2,4,5], 4, 7),
 ('sigma (r,s)^3 d7 [triple along (01)-fixed surface]', Msig, [1,4], 3, 7),
 ('sigma (p,q,b)^3 d7', Msig, [0,3,5], 3, 7),
 ('sigma (r,s,a)^3 d7 [elliptic cone along curve X∩{X0=X1,Y0=Y1,X2=0}]', Msig, [1,4,2], 3, 7),
 ('sigma (r,s,b)^3 d7', Msig, [1,4,5], 3, 7),
 ('orig (X0,X1,Y0,Y1)^4 d7 [quadruple pts on (X2,Y2)-line]', None, [0,1,3,4], 4, 7),
 ('sigma (r,s)^2 d7 [double along (01)-fixed surface]', Msig, [1,4], 2, 7),
 ('orig (X2,Y2)^2 d7 [double along X∩{X2=Y2=0}]', None, [2,5], 2, 7),
 ('orig (X0,X1,X2,Y0)^4 d7', None, [0,1,2,3], 4, 7),
]
allfam = {}
for label, M, idx, k, d in specs:
    reps = family(idx, k, d, verbose=False, Mx=M, My=M)
    allfam[label] = [to_str(P) for P in reps]
    n = len(reps)
    print(f'== {label}: dim {n}', flush=True)
    if n == 0: continue
    if n <= 3:
        for i, P in enumerate(reps): run(f'{label} basis#{i}', P, label)
        if n >= 2:
            run(f'{label} sum of basis', add(reps[0], reps[1]), label)
            run(f'{label} difference of basis', add(reps[0], scale(reps[1], Fr(-1))), label)
    else:
        inv = invariant_part(reps); ti = t01_invariant(reps); ta = t01_anti(reps)
        print(f'   S3-inv dim {len(inv)}, (01)-inv dim {len(ti)}, (01)-anti dim {len(ta)}', flush=True)
        if inv: run(f'{label} S3-invariant random combo', rnd(inv), label)
        if ti: run(f'{label} (01)-invariant random combo', rnd(ti), label)
        if ta: run(f'{label} (01)-anti random combo', rnd(ta), label)
        run(f'{label} generic random combo', rnd(reps), label)
json.dump(allfam, open('/tmp/hodge/results/diagonal-permutation/families_d7_adapted.json', 'w'), indent=1)
