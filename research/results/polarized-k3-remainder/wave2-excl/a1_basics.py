"""(a) re-verification 1: exact identities, Hilbert functions over Q, eta != 0, UV-frame symmetry group."""
import sys, time, json
sys.path.insert(0, '.')
from k3lib import *
from logutil import log
from hodge import symmetry
t0 = time.time()
F, A, q = load_xy()
Fh, Ah, _ = load_uv()
f = ternary_f()
# 1. F = (f(U)+f(V))/2 exactly, A-hat = 2^-6 prod (U_i^2 - V_i^2)
fU = {(e[0], e[1], e[2], 0, 0, 0): c for e, c in f.items()}
fV = {(0, 0, 0, e[0], e[1], e[2]): c for e, c in f.items()}
chk1 = sub(Fh, scale(add(fU, fV), Fr(1, 2))) == {}
prod = {(0,) * NV: Fr(1, 64)}
for i in range(3):
    eu = [0] * NV; eu[i] = 2; ev = [0] * NV; ev[3 + i] = 2
    prod = mul(prod, {tuple(eu): Fr(1), tuple(ev): Fr(-1)})
chk2 = sub(Ah, prod) == {}
# round trip
chk3 = sub(uv_to_xy(Fh), F) == {} and sub(uv_to_xy(Ah), A) == {}
# parities
Fm = {e: c * (-1) ** (e[3] + e[4] + e[5]) for e, c in F.items()}
Am = {e: c * (-1) ** (e[3] + e[4] + e[5]) for e, c in A.items()}
chk4 = sub(Fm, F) == {} and add(Am, A) == {}
print('F=(f(U)+f(V))/2:', chk1, '| Ahat = 2^-6 prod(U_i^2-V_i^2):', chk2, '| round trip:', chk3, '| parities:', chk4)
print('f =', {e: str(c) for e, c in sorted(f.items())})
# 2. F mod 5 is Fermat
F5 = {e: Fr(c.numerator * pow(c.denominator, -1, 5) % 5) for e, c in F.items() if (c.numerator * pow(c.denominator, -1, 5)) % 5}
fermat = {tuple(6 if j == i else 0 for j in range(NV)): Fr(1) for i in range(NV)}
chk5 = F5 == fermat
print('F mod 5 == Fermat sextic:', chk5)
# 3. Hilbert function of R_f over Q, eta = [z0^2 z1^2 z2^2] != 0
R = TernaryRingQ(f, maxdeg=14)
hilb = [R.dim(k) for k in range(15)]
print('Hilb(R_f) over Q:', hilb)
eta = R.nf({(2, 2, 2): Fr(1)}, 6)
print('eta normal form (nonzero entries):', len(eta), 'eta != 0:', bool(eta))
# 4. tensor ring dims (should be the complete-intersection Hilbert series (1+t+..+t^4)^6)
T = TensorRingQ(R, maxdeg=14)
dims = [T.dim(d) for d in range(15)]
import sympy
t = sympy.symbols('t')
ci = sympy.Poly(sympy.expand(sum(t ** i for i in range(5)) ** 6), t).all_coeffs()[::-1]
print('dim R_F,d over Q (tensor):', dims)
print('expected (1+..+t^4)^6:', [int(c) for c in ci[:15]])
chk6 = dims == [int(c) for c in ci[:15]]
# 5. UV-frame diagonal symmetry group of Fhat and characters of Ahat's components
G = symmetry.diagonal_group(Fh)
print('UV-frame projective diagonal group order:', G.order, 'elements:', G.elements)
Gxy = symmetry.diagonal_group(F)
print('XY-frame projective diagonal group order:', Gxy.order)
comps = split_vdeg_mod(Ah)
chars = {}
for r, P in comps.items():
    chars[r] = [symmetry.character(P, w, G.order) for w in G.elements]
    rc = [symmetry.residue_character(Fh, P, 1, w, G.order) for w in G.elements]
    print(f'Ahat component V-degree={r} mod 6: {len(P)} monomials, characters {chars[r]}, residue chars {rc}')
log(exp='a1_basics', checks=dict(F_eq_half_fU_fV=chk1, Ahat_product=chk2, roundtrip=chk3, parities=chk4, F_mod5_Fermat=chk5,
    tensor_dims_match_CI=chk6), hilb_Rf_Q=hilb, eta_nonzero=bool(eta), dims_RF_Q=dims,
    uv_group_order=G.order, uv_group_elements=[list(w) for w in G.elements], xy_group_order=Gxy.order,
    ahat_components={str(k): len(v) for k, v in comps.items()}, seconds=round(time.time() - t0, 1), rigor='exact over Q')
