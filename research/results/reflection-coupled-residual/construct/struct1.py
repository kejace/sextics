# Structural checks for the CONSTRUCTION lens (sympy, exact).
# 1. P in I_X * I_Y (every coupling monomial has >=2 distinct X and >=2 distinct Y variables) -> six Fermat surfaces.
# 2. Exchange-cover scan: for every coordinate triple (l1,l2 | l3) test the symmetry criterion h_{c,m} = 0 for c != m.
# 3. Galois closure data: a, b, discriminant; squarefreeness/irreducibility checks.
import sympy as sp, itertools, json, time
X0,X1,X2,Y0,Y1,Y2 = Z = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
c = sp.Rational(5,7)
P = X0**3*X1*Y0*Y1 + X0*X1*Y0**3*Y1 + X0*X1*Y0**2*Y2**2 + 2*X0*X1*Y1**2*Y2**2 + 3*X0*X2**3*Y0*Y1 + 2*X1**3*X2*Y0*Y1 + X1*X2*Y0*Y1**3
F = sum(z**6 for z in Z) + c*P
A = X0**3*Y0**3
out = {}
# 1.
Pp = sp.Poly(P, *Z)
ok = True
for mon, coef in Pp.terms():
    nx = sum(1 for e in mon[:3] if e > 0); ny = sum(1 for e in mon[3:] if e > 0)
    if nx < 2 or ny < 2: ok = False
out['P_in_IX_IY'] = ok
# the six Fermat surfaces: F restricted to Z_i = Z_j = 0 (i,j same block)
surf = {}
for blk in [(0,1,2),(3,4,5)]:
    for i,j in itertools.combinations(blk,2):
        Fr = F.subs({Z[i]:0, Z[j]:0})
        surf[f'{Z[i]}={Z[j]}=0'] = str(sp.expand(Fr))
out['coordinate_plane_sections'] = surf
# mixed pairs (one X one Y): is the section still Fermat?
mixed = {}
for i in range(3):
    for j in range(3,6):
        Fr = sp.expand(F.subs({Z[i]:0, Z[j]:0}))
        mixed[f'{Z[i]}={Z[j]}=0'] = 'Fermat' if Fr == sum(z**6 for z in Z if z not in (Z[i],Z[j])) else str(Fr)
out['mixed_plane_sections'] = mixed
# 2. exchange-cover scan.  Block (l1,l2 | l3), complement Zp.  F = sum_{c,m} l3^c h_{c,m}(l1,l2; Zp), h_{c,m} of degree m in (l1,l2).
# Criterion (proved in the note): an exchange T <-> l3 over the base curve g(alpha,beta)=1 (g = F|_{span(l1,l2)})
# exists only if h_{c,m} = 0 for all c != m with (c,m) != (6,0),(0,6).  Report the offending (c,m) list.
scan = []
for l3 in range(6):
    for l1,l2 in itertools.combinations([k for k in range(6) if k != l3], 2):
        Fp = sp.Poly(F, *Z)
        bad = {}
        for mon, coef in Fp.terms():
            cdeg = mon[l3]; mdeg = mon[l1] + mon[l2]
            if cdeg != mdeg and (cdeg, mdeg) not in [(6,0),(0,6)]:
                bad.setdefault(f'({cdeg},{mdeg})', []).append(str(sp.Mul(*[Z[k]**mon[k] for k in range(6)])))
        scan.append({'l1l2': [str(Z[l1]), str(Z[l2])], 'l3': str(Z[l3]), 'n_bad_types': len(bad), 'bad': bad})
out['exchange_scan'] = scan
out['exchange_scan_any_pass'] = any(s['n_bad_types'] == 0 for s in scan)
# 3. Galois closure of X -> P^4_x (forget Y2). F = u^3 + a u + b, u = Y2^2.
u = sp.Symbol('u')
H = sp.expand(F.subs(Y2, sp.sqrt(u)))
Hp = sp.Poly(H, u)
out['H_coeffs_in_u'] = [str(sp.factor(cf)) for cf in Hp.all_coeffs()]
a = Hp.coeff_monomial(u); b = Hp.coeff_monomial(1)
Delta = sp.expand(-4*a**3 - 27*b**2)
x5 = (X0,X1,X2,Y0,Y1)
out['a_factor'] = str(sp.factor(a)); out['b_nterms'] = len(sp.Poly(b,*x5).terms())
# squarefree / irreducible tests over Q
fb = sp.factor_list(b); fD = sp.factor_list(Delta)
out['b_factors'] = [(str(f)[:60], e) for f,e in fb[1]]
out['Delta_factors_degrees'] = [(sp.Poly(f,*x5).total_degree(), e) for f,e in fD[1]]
out['Delta_is_square'] = all(e % 2 == 0 for f,e in fD[1])
out['b_is_square'] = all(e % 2 == 0 for f,e in fb[1])
# Galois closure Zhat: (2,4,6) CI in P^7 (x, y1,y2,y3): y1^2+y2^2+y3^2 = 0, sum y_i^2 y_j^2 = a, y1^2 y2^2 y3^2 = -b.
out['Zhat'] = 'CI of degrees (2,4,6) in P^7, K = O(2+4+6-8) = O(4); singular along {y=0, a=b=0} (A1 x C^2 type, canonical)'
json.dump(out, open('/tmp/hodge/results/reflection-coupled-residual/construct/struct1.json','w'), indent=1)
print(json.dumps({k:v for k,v in out.items() if k not in ('exchange_scan','coordinate_plane_sections')}, indent=1))
print('exchange scan: triples with zero offending types:', [s for s in scan if s['n_bad_types']==0])
print('min offending types:', min(s['n_bad_types'] for s in scan))
for s in scan[:6]: print(s['l1l2'], s['l3'], s['n_bad_types'], list(s['bad'].keys()))
