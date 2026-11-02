"""Cones over plane curves with vertex plane meeting X in 6 concurrent lines (vertex point on S_H), mod p."""
from h import *
import random
random.seed(5)
p = core.PRIMES[0]
t = 388111789      # root of f(t,1) mod p  -> vertex [t:1:0:0:0:0] on S_H
s = 1961481078     # root of g(s,1) mod p  -> vertex [0:0:s:1:0:0] on S_H
Fp = core.reduce_mod(F, p)
l = {(1,0,0,0,0,0): 1, (0,1,0,0,0,0): (-t) % p}      # X0 - t X1
m = {(0,0,0,1,0,0): 1, (0,0,0,0,1,0): (-s) % p}      # Y0 - s Y1
lp = {(1,0,0,0,0,0): 1, (0,1,0,0,0,0): (-(t+1)) % p} # control: vertex off S_H
Y0 = P('Y0'); Y1 = P('Y1'); X0 = P('X0'); X1 = P('X1')
g = core.reduce_mod(P('Y0**6 + 5*Y0**5*Y1/7 + 15*Y0*Y1**5/7 + Y1**6'), p)
f = core.reduce_mod(P('X0**6 + 5*X0**5*X1/7 + 10*X0*X1**5/7 + X1**6'), p)
def pw(Q, k): return core.power(Q, k, p)
def addp(*Qs):
    out = {}
    for Q in Qs: out = core.add(out, Q, p)
    return out
def sc(Q, c): return core.scale(Q, c % p, p)
# sanity: f(t,1)=0 mod p
assert core.evaluate(f, (t,1,0,0,0,0), p) == 0 and core.evaluate(g, (0,0,0,s,1,0), p) == 0
tests = []
tests.append((addp(pw(l,6), pw(Y0,6), pw(Y1,6)), 'cone3-vertex-on-SH', 'l^6+Y0^6+Y1^6, l=X0-tX1, f(t,1)=0'))
tests.append((addp(pw(l,6), g), 'cone3-vertex-on-SH', 'l^6+g(Y0,Y1), vertex on S_H'))
tests.append((addp(pw(l,6), sc(g, -1)), 'cone3-vertex-on-SH', 'l^6-g(Y0,Y1), vertex on S_H'))
tests.append((addp(pw(lp,6), pw(Y0,6), pw(Y1,6)), 'cone3-vertex-off-SH(control)', "l'^6+Y0^6+Y1^6, l'=X0-(t+1)X1 (f(t+1,1)!=0)"))
tests.append((addp(pw(X0,6), pw(X1,6), pw(m,6)), 'cone3-vertex-on-SH', 'X0^6+X1^6+m^6, m=Y0-sY1, g(s,1)=0'))
tests.append((addp(f, pw(m,6)), 'cone3-vertex-on-SH', 'f(X0,X1)+m^6, vertex on S_H'))
# random even sextic in (l, Y0, Y1)
for k in range(2):
    D = {}
    for a in range(7):
        for b in range(7 - a):
            c = 6 - a - b
            if (b + c) % 2: continue
            D = core.add(D, sc(core.mul(core.mul(pw(l, a), pw(Y0, b), p), pw(Y1, c), p), random.randint(1, 50)), p)
    tests.append((D, 'cone3-vertex-on-SH-random', f'random even sextic in (l,Y0,Y1) #{k}, vertex on S_H'))
# cubic and quartic cones with vertex on S_H
tests.append((addp(pw(Y0,3), pw(Y1,3), core.mul(pw(l,2), Y0, p)), 'cone3-vertex-on-SH', 'Y0^3+Y1^3+l^2 Y0 (Y-odd cubic), vertex on S_H'))
tests.append((addp(pw(l,4), pw(Y0,4), pw(Y1,4)), 'cone3-vertex-on-SH', 'l^4+Y0^4+Y1^4, vertex on S_H'))
tests.append((addp(pw(X0,4), pw(X1,4), pw(m,4)), 'cone3-vertex-on-SH', 'X0^4+X1^4+m^4, vertex on S_H'))
order = [8, 9, 0, 1, 3, 4, 5, 2, 6, 7, 10]
for i in order:
    D, fam, note = tests[i]
    r = run_modp(D, p, fam, expr=note, note=note)
    print({k: r[k] for k in ('D','d','consistent','rank','nullity','rows','cols','seconds')}, flush=True)
