import json
R = json.load(open('analysis.json'))
def z(u):
    a, b = u
    if b == 0: return f"{a}"
    return f"{a}{'+' if b>0 else '-'}{abs(b)}z"
print("## split primes")
print("| p | #T_B(F_p) | #E'(F_p) | pi_E' | b1 | b2 | b3 | b3 cross-check F_p^3 | Newton at pp | L_P factorisation over Q | ratio root of unity (exact/numeric) | j=0 roots in L_+ |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
simple = []; cmroots = []; classes = {}
for r in R:
    if r['class'] != 'split': continue
    if r.get('b3') is None and r.get('b3_from_Fp3') is None:
        print(f"| {r['p']} | {r['#TB(Fp)']} | {r['#E(Fp)']} | {z(r['piE'])} | 0 | 0 | undetermined | | | | | |"); continue
    b3 = r['b3'] if r['b3'] is not None else r['b3_from_Fp3']
    cc = r.get('b3_crosscheck', '')
    fac = ' * '.join(f"({f})^{e}" if e > 1 else f"({f})" for f, e in r['factor_Q'])
    ru = f"{r['rootunity_hits']}/{r['rootunity_numeric']}"
    print(f"| {r['p']} | {r['#TB(Fp)']} | {r['#E(Fp)']} | {z(r['piE'])} | {z(r['b1'])} | {z(r['b2'])} | {z(b3)} | {cc} | {r['newton_pp(L+P)']} | {fac} | {ru} | {[z(g) for g in r['cm_elliptic_roots_in_L+']]} |")
    if len(r['factor_Q']) == 1 and r['factor_Q'][0][1] == 1 and not r['rootunity_hits'] and not r['rootunity_numeric']:
        simple.append(r['p'])
        dc = r['disc_class']; classes[r['p']] = (tuple(sorted((e, tuple(g)) for e, g in dc['odd_primes'])), dc['unit_sign'])
    if r['cm_elliptic_roots_in_L+']: cmroots.append(r['p'])
print()
print("primes where P mod p is geometrically simple (L_P irreducible sextic, no root-of-unity ratio):", simple)
print("primes where L_+ has a root in Z[zeta] of norm p (j=0 elliptic factor of the reduction):", cmroots)
print("relative discriminant classes of L_+ over Q(zeta) (primes with odd exponent, unit sign):")
for p in sorted(classes): print(p, classes[p])
print()
print("## inert primes")
print("| p | #T_B(F_p^2) | #E'(F_p^2) | E' supersingular | t2 real | M(S) | c/p | second factor j=0-type | second factor ordinary | F_p^4 check |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in R:
    if r['class'] != 'inert': continue
    print(f"| {r['p']} | {r['#TB(Fp2)']} | {r['#E(Fp2)']} | {r['E_supersingular_check']} | {r['t2_real']} | {r['M(S)']} | {r['c_over_p']} | {r['second_factor_j0_type']} | {r['second_factor_ordinary']} | {r.get('Fp4_check','')} |")
