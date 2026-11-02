"""Analyse the zero set of an L-relaxation after square elimination: residual homogeneous quadratic ideal on the
surviving parameters -> Singular (dimension, degree, minimal associated primes when feasible) over F_p.
Also prints the surviving parameters as (i, monomial) so components can be read off as denominator families.
usage: survivors.py <relax json prefix> [--primes] [--charts]"""
import sys, json, subprocess, os, time
import numpy as np
prefix = sys.argv[1]
J = json.load(open(prefix + '.json'))
Z = np.load(prefix + '.npz')
C = Z['C']; pairs = [tuple(x) for x in Z['pairs']]
p = J['p']; pars = J['parameter_basis']
names = ['U0', 'U1', 'U2', 'V0', 'V1', 'V2']
def mon_str(m): return '*'.join(f'{names[i]}^{e}' if e > 1 else names[i] for i, e in enumerate(m) if e) or '1'
killed = set(J.get('square_elimination_mod_p', {}).get('killed', []))
surv = [i for i in range(len(pars)) if i not in killed]
print('params', len(pars), 'killed', len(killed), 'survivors', len(surv))
for i in surv: print('  c%d = coefficient of %s in a_%d (P_%d * %s)' % (i, mon_str(pars[i][1]), pars[i][0], pars[i][0], mon_str(pars[i][1])))
keep = [j for j, (a, b) in enumerate(pairs) if a in surv and b in surv]
sub = np.ascontiguousarray(C[:, keep]) % p
from flint import nmod_mat
M = nmod_mat(sub.shape[0], sub.shape[1], sub.ravel().tolist(), p)
rr, rank = M.rref()
print('residual constraint rank', rank, 'on', len(keep), 'quadratic monomials in', len(surv), 'variables')
polys = []
for i in range(rank):
    terms = []
    for jj in range(len(keep)):
        v = int(rr[i, jj])
        if v:
            a, b = pairs[keep[jj]]
            terms.append(f'{v}*c{a}*c{b}' if a != b else f'{v}*c{a}^2')
    polys.append('+'.join(terms))
vars_ = ','.join(f'c{i}' for i in surv)
sing = [f'ring r = {p}, ({vars_}), dp;', 'option(redSB);', 'LIB "primdec.lib";', 'ideal I = ' + ',\n'.join(polys) + ';',
        'ideal G = std(I);', 'print("DIM=" + string(dim(G)) + " DEG=" + string(mult(G)) + " SIZE=" + string(size(G)));']
if '--primes' in sys.argv:
    sing += ['list L = minAssGTZ(I);', 'print("NCOMP=" + string(size(L)));', 'int k; for (k=1; k<=size(L); k++) { print("COMPONENT " + string(k) + " dim " + string(dim(std(L[k])))); print(L[k]); }']
sing.append('quit;')
sf = prefix + '-residual.sing'
open(sf, 'w').write('\n'.join(sing) + '\n')
t0 = time.time()
out = subprocess.run(['Singular', '-q', sf], capture_output=True, text=True, timeout=3600)
print(out.stdout[-6000:]); print(out.stderr[-500:])
print('singular seconds', round(time.time() - t0, 1))
open(prefix + '-residual.log', 'w').write(out.stdout + out.stderr)
