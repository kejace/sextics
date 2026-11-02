import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *
iid = 'coupled-pair-rational-d4-strict'
F, A, q, rec = core.load_instance(iid)
d, D, N, Q = core.load_witness(iid)
print('A =', core.to_str(A), ' pA =', parity_of(A), ' pD =', parity_of(D), ' Q terms', len(Q))
N2 = normalize_lift(F, A, D, N)
for i in range(6):
    print(i, 'N parity', parity_of(N[i]), 'N2 parity', parity_of(N2[i]), 'deg', core.degree(N2[i]), 'terms', len(N2[i]))
prob = E2Float(F, A, d, parity_of(D), verbose=True)
print('pN', prob.pN, 'SdivN par', prob.SdivN.par)
Dv = prob.SD.vec(D); Nv = prob.Nvec(N2)
# exact pieces
divN = {}
for i in range(6): divN = add(divN, pderiv(N2[i], i))
print('exact divN terms', len(divN), 'parity', parity_of(divN))
divNf = np.zeros(prob.SdivN.n)
for i in range(6): divNf += prob.dN[i] @ Nv[prob.Noff[i]:prob.Noff[i+1]]
print('float divN vs exact:', np.linalg.norm(divNf - prob.SdivN.vec(divN)))
P1 = mul(D, divN)
P1f = prob.T_DdivN.prod(Dv, divNf)
print('D*divN float vs exact:', np.linalg.norm(P1f - prob.Srow.vec(P1)), np.linalg.norm(P1f))
for i in range(6):
    Di = pderiv(D, i); Dif = prob.dD[i] @ Dv
    print(i, 'D_i float vs exact', np.linalg.norm(Dif - prob.SDi[i].vec(Di)) if Di else np.linalg.norm(Dif),
          'N_i D_i:', np.linalg.norm(prob.T_NiDi[i].prod(Nv[prob.Noff[i]:prob.Noff[i+1]], Dif) - prob.Srow.vec(mul(N2[i], Di))))
B = prob.Bvec(Dv, Nv); Bex = B_of(D, N2)
print('|Bfloat|', np.linalg.norm(B), '|Bexact|', len(Bex))
