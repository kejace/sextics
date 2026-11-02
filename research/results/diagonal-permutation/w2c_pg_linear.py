"""J3: p_g of the quotient of X cap P^3 by G=<s,tau> for every G-stable P^3 (sum of eigenspaces), s=(01), tau=Y->-Y.
Eigen-coordinates: p=X0+X1, a=X2 (char ++), q=Y0+Y1, b=Y2 (+-), r=X0-X1 (-+), t=Y0-Y1 (--).  Sextic surface S in P^3,
K_S = O(2): 2-forms = quadrics G with character chi_G = det(G on P^3)^{-1} (G-invariant forms = p_g of S/G)."""
import itertools, json
chars = {'p':(1,1),'a':(1,1),'q':(1,-1),'b':(1,-1),'r':(-1,1),'t':(-1,-1)}
# types of stable 4-dim subspaces: choose a 'generic' member of each eigenspace when a 1-dim piece is taken
types = {
 'A: E++ + E+-  (fixed P^3 of s)': ['p','a','q','b'],
 'B: E++ + E-+ + E--': ['p','a','r','t'],
 'C: E++ + l + E-+ (l in E+-)': ['p','a','q','r'],
 'D: E++ + l + E--': ['p','a','q','t'],
 'E: l\' + E+- + E-+': ['p','q','b','r'],
 'F: l\' + E+- + E--': ['p','q','b','t'],
 'G: l\' + l + E-+ + E-- (2-param family)': ['p','q','r','t'],
}
out = {}
for name, coords in types.items():
    det = (1,1)
    for c in coords: det = (det[0]*chars[c][0], det[1]*chars[c][1])
    inv = []
    for i,j in itertools.combinations_with_replacement(coords, 2):
        ch = (chars[i][0]*chars[j][0], chars[i][1]*chars[j][1])
        if ch == det: inv.append(i+j)   # chi_G * det = 1  <=> chi_G = det (characters are +-1)
    out[name] = (len(inv), inv)
    print(name, '-> p_g(S/G) =', len(inv), inv)
with open('wave2-construct.log.jsonl','a') as f:
    f.write(json.dumps({"exp":"J3-pg-linear-sections","family":"<s,tau>-stable P^3 sections of X, quotient p_g","result":{k:v[0] for k,v in out.items()},"job":"local","note":"p_g(S/G)=#quadrics G on P^3 with chi_G = det; min 1 (types B-F, 1-parameter), 4 for the only 2-parameter family (G)"})+"\n")
