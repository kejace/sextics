import sys, json
sys.path.insert(0, '/tmp/hodge/results/eight-plane-configuration')
from lib import *
import lib
_log = lib.log
def log2(e): e['engine'] = 2; _log(e)
lib.log = log2
K = [core.parse(s) for s in json.load(open(f'{WD}/K_decomp.json'))['K']]
H = [core.parse(s) for s in json.load(open(f'{WD}/K_decomp.json'))['H']]
Q = [Q0, Q1, Q2]
# verify decomposition F = sum Q_i K_i
chk = {}
for i in range(3): chk = core.add(chk, core.mul(Q[i], K[i]))
assert core.sub(chk, F) == {}, 'decomp mismatch'
print('F = sum Q_i K_i verified exactly')
diag = [P(f'X{i}**4 - X{i}**2*Y{i}**2 + Y{i}**4') for i in range(3)]
ferm = P('X0**6+X1**6+X2**6+Y0**6+Y1**6+Y2**6')
cands = []
for i in range(3):
    cands += [(f'Q{i}', 'Q-products', Q[i]), (f'Q{i}**2', 'Q-products', core.mul(Q[i],Q[i])), (f'Q{i}**3','Q-products', core.power(Q[i],3)),
              (f'K{i}', 'K-quartics', K[i]), (f'H{i}', 'H-quartics', H[i]), (f'diag{i}', 'diag-quartic', diag[i]),
              (f'Q{i}*K{i}', 'block-terms', core.mul(Q[i],K[i])), (f'F - Q{i}*K{i}', 'partial-sums', core.sub(F, core.mul(Q[i],K[i]))),
              (f'X{i}**6+Y{i}**6', 'pure-power', P(f'X{i}**6+Y{i}**6')), (f'Q{i}*diag{i}', 'pure-power', core.mul(Q[i],diag[i]))]
for i in range(3):
    for j in range(i+1,3):
        cands += [(f'Q{i}*Q{j}', 'Q-products', core.mul(Q[i],Q[j])), (f'Q{i}*K{j}+Q{j}*K{i}','cross-terms', core.add(core.mul(Q[i],K[j]), core.mul(Q[j],K[i]))),
                  (f'Q{i}**2*Q{j}', 'Q-products', core.mul(core.mul(Q[i],Q[i]),Q[j]))]
cands += [('Q0*Q1*Q2', 'Q-products', core.mul(core.mul(Q0,Q1),Q2)), ('Q0+Q1+Q2', 'Q-linear', core.add(core.add(Q0,Q1),Q2)),
          ('(Q0+Q1+Q2)**2','Q-linear', core.power(core.add(core.add(Q0,Q1),Q2),2)), ('(Q0+Q1+Q2)**3','Q-linear', core.power(core.add(core.add(Q0,Q1),Q2),3)),
          ('Fermat', 'pure-power', ferm), ('F - Fermat', 'pure-power', core.sub(F, ferm)),
          ('X0**3+X1**3+X2**3', 'cone-plane-vertex', P('X0**3+X1**3+X2**3')), ('Y0**3+Y1**3+Y2**3', 'cone-plane-vertex', P('Y0**3+Y1**3+Y2**3')),
          
          ('X0**4+X1**4+Y0**4+Y1**4', 'cone-line-vertex', P('X0**4+X1**4+Y0**4+Y1**4')), ('Q0**2+Q1**2', 'cone-line-vertex', core.add(core.mul(Q0,Q0),core.mul(Q1,Q1))),
          ('X0**4+Y0**4+X1**2*Y1**2', 'cone-line-vertex', P('X0**4+Y0**4+X1**2*Y1**2')),
          ('K0+K1+K2', 'K-quartics', core.add(core.add(K[0],K[1]),K[2])), ('H0+H1+H2', 'H-quartics', core.add(core.add(H[0],H[1]),H[2])),
          ('sum Q_i*diag_i', 'pure-power', core.add(core.add(core.mul(Q0,diag[0]),core.mul(Q1,diag[1])),core.mul(Q2,diag[2]))),
          ('sum Q_i*H_i', 'H-quartics', core.add(core.add(core.mul(Q0,H[0]),core.mul(Q1,H[1])),core.mul(Q2,H[2]))),
          ('X0*Y0', 'odd-low', P('X0*Y0')), ('X0*Y0*Q0','odd-low', core.mul(P('X0*Y0'),Q0)), ('X0**3*Y0**3','odd-low', P('X0**3*Y0**3')),
          ('X0*Y0*Q1*Q2','odd-low', core.mul(P('X0*Y0'),core.mul(Q1,Q2))), ('X0*Y0*K0','odd-low', core.mul(P('X0*Y0'),K[0]))]
start = [i for i,c in enumerate(cands) if c[0]=='X0**4+X1**4+Y0**4+Y1**4'][0]
for name, fam, D in cands[start:]:
    e = test_D(D, name, fam, note='geometric family; E1 only (kernel theorem predicts failure unless D in J)', e1_only=True)
    print(f"{name:28s} d={e['d']} e1={e['e1']} {e['e1_info']} {e['seconds']}s", flush=True)
