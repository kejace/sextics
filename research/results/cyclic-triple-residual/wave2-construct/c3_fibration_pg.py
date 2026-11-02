# C3: surface-fibration companions of X/K, K inside the stabiliser <s> of the class (s=(1,5,5,2,2,0), s^3 = tau).
# Linear fibrations X -> P^2 (base = 3 coordinates B, fibre S_b = sextic surface in P^3 with coordinates lambda + 3 others).
# K_b = {s^k : the three base coordinates have equal s^k-weight} acts on the fibre; p_g(S_b/K_b) = number of K_b-invariant
# residues m*Omega_3/F_b with m quadratic in the fibre coordinates (Griffiths), weight(residue) = wt(m) + wt(vol form).
import itertools, json
names = ['X0','X1','X2','Y0','Y1','Y2']
s = [1,5,5,2,2,0]
out = []
for B in itertools.combinations(range(6), 3):
    fib = [i for i in range(6) if i not in B]
    # fibre coordinates: lambda (weight = common base weight w_B if equal) and the 3 fibre coordinates
    Kb = []
    for k in range(6):
        wts = {(s[i]*k) % 6 for i in B}
        if len(wts) == 1: Kb.append(k)
    # weights of fibre coords under s^k: lambda has the base weight
    res = []
    quad = list(itertools.combinations_with_replacement(range(4), 2))  # fibre coords 0=lambda,1..3
    def wt(k, c):
        return (s[B[0]]*k) % 6 if c == 0 else (s[fib[c-1]]*k) % 6
    pg = 0; inv = []
    for m in quad:
        ok = True
        for k in Kb:
            vol = sum(wt(k, c) for c in range(4)) % 6
            if (wt(k, m[0]) + wt(k, m[1]) + vol) % 6 != 0: ok = False; break
        if ok: pg += 1; inv.append(m)
    out.append({'base': [names[i] for i in B], 'K_b_exponents': Kb, 'order_Kb': len(Kb), 'pg_fibre_quotient': pg})
for o in sorted(out, key=lambda o: (o['pg_fibre_quotient'], -o['order_Kb'])):
    print(o)
json.dump(out, open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c3_fibration_pg.json','w'), indent=1)
print('minimum p_g over all equivariant linear P^2-fibrations:', min(o['pg_fibre_quotient'] for o in out))
