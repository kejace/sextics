"""Fermat calibration: which denominators serve X0^3Y0^3 on the Fermat sextic?"""
import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
import ctr
from ctr import *
F0 = P('X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6')
A0 = A
out = []
def run(Dexpr, fam, note=''):
    r = test_D(Dexpr, fam, note='FERMAT-CALIBRATION ' + note, F_=F0, A_=A0, q_=1, tag='fermat')
    print(fmt(r), flush=True); out.append(r); return r
run('X0**6 + Y0**6', 'fermat:cone-hyperplanes', 'six cone hyperplanes X0=lambda Y0, lambda^6=-1')
run('X0**3 + Y0**3', 'fermat:half', 'three hyperplanes, not cones')
run('X0**2 + Y0**2', 'fermat:quadric')
run('X0**6 - Y0**6', 'fermat:cone-hyperplanes-wrong-sign')
run('X0**3 + X1**3 + X2**3', 'fermat:cone3-X', 'cubic cone over elliptic curve, vertex P^2_Y')
run('Y0**3 + Y1**3 + Y2**3', 'fermat:cone3-Y', 'cubic cone, vertex P^2_X')
run('X0**3 + X1**3 + X0*Y0**2', 'fermat:cone3-mixed')
run('Y0**3 + Y1**3 + X0**2*Y0', 'fermat:cone3-mixed-odd')
run('X0**4 + X1**4 + X2**4 + Y0**4', 'fermat:cone4-quartic', 'K3 base, vertex line Y1Y2')
run('X0**6 + X1**6 + X2**6 + Y0**6', 'fermat:cone4-sextic', 'sextic base p_g=10, vertex line Y1Y2')
run('X1**6 + X2**6 + Y1**6 + Y2**6', 'fermat:sum4-complement', 'complement of {X0,Y0}: join structure')
run('X0**6 + X1**6 + X2**6', 'fermat:sum3-X', 'join of two plane sextics')
run('X0**6 + Y0**6 + Y1**6', 'fermat:sum3')
run('X1**6 + X2**6', 'fermat:sum2')
run('X0**6 + 2*Y0**6', 'fermat:sum2-coef2')
run('X0**4 + Y0**4', 'fermat:quartic-binary')
run('X0**5 + Y0**5', 'fermat:quintic-binary')
json.dump(out, open('batch_fermat.json', 'w'), indent=1)
