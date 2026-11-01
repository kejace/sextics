import sys, time, json, itertools
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual')
from j6 import *
def G_diag(m):
    return [core.scale(Zlin[i], Fraction(m[i])) if m[i] else {} for i in range(6)]
# D_* and close relatives
tests = [
  ((0,0,0,1,-1,1), 'D* = Y0F_Y0 - Y1F_Y1 + Y2F_Y2 (quadruple curve along C_X)'),
  ((0,0,0,1,1,1),  'E_Y(F)'),
  ((1,1,1,0,0,0),  'E_X(F)'),
  ((0,0,0,1,-1,-1), 'Y0F_Y0 - Y1F_Y1 - Y2F_Y2'),
  ((0,0,0,-1,1,1), '-Y0F_Y0 + Y1F_Y1 + Y2F_Y2'),
  ((0,0,0,1,1,-1), 'Y0F_Y0 + Y1F_Y1 - Y2F_Y2'),
  ((0,0,1,0,0,0),  'X2F_X2 = 6 X2^3 U'),
  ((1,0,0,0,0,0),  'X0F_X0'),
  ((0,1,0,0,0,0),  'X1F_X1'),
  ((0,0,0,1,0,0),  'Y0F_Y0'),
  ((0,0,0,0,1,0),  'Y1F_Y1'),
  ((0,0,0,0,0,1),  'Y2F_Y2'),
  ((1,1,1,1,1,1),  '6F (illegal control)'),
]
for m, name in tests:
    ent, r = test_G(G_diag(m), f'diag{m}: {name}', family='J6-diag', note='weighted Euler derivative sum m_i Z_i F_i')
