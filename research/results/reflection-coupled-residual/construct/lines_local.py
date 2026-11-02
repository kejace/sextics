# Local dimension of F_1(X) at representative Fermat lines (tangent dims 4,3,2,1), via Singular local ordering (ds).
import sys, itertools, subprocess, json, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
from hodge import core
p = 65521
F, A, q, rec = core.load_instance('reflection-coupled-residual')
def modp(fr): return (fr.numerator % p) * pow(fr.denominator % p, p-2, p) % p
roots = [z for z in range(1, p) if pow(z, 6, p) == p-1]
Z = ['X0','X1','X2','Y0','Y1','Y2']
def line(blk_pair, pairing_idx, z1, z2):
    rest = [j for j in range(6) if j not in blk_pair]
    pairing = [((rest[0],rest[1]),(rest[2],rest[3])), ((rest[0],rest[2]),(rest[1],rest[3])), ((rest[0],rest[3]),(rest[1],rest[2]))][pairing_idx]
    (k,l),(m,n) = pairing
    pv = [0]*6; qv = [0]*6; pv[k] = z1; pv[l] = 1; qv[m] = z2; qv[n] = 1
    return pv, qv, (k,l,m,n)
def singular_local_dim(pv, qv, rs):
    r, s = rs  # chart indices with p_r=1,p_s=0,q_r=0,q_s=1  (we pick r=l, s=n)
    free = [j for j in range(6) if j not in (r, s)]
    # coordinates: Z_j = (pv_j + a_j) + t (qv_j + b_j) for j free; Z_r = 1, Z_s = t
    Fstr = ' + '.join(f'{modp(c)}*' + '*'.join(f'{Z[j]}^{e[j]}' for j in range(6) if e[j]) for e, c in F.items())
    subs = []
    for j in range(6):
        if j == r: subs.append(f'{Z[j]}, 1')
        elif j == s: subs.append(f'{Z[j]}, t')
        else: subs.append(f'{Z[j]}, ({pv[j]} + a{j}) + t*({qv[j]} + b{j})')
    avars = ','.join(f'a{j}' for j in free); bvars = ','.join(f'b{j}' for j in free)
    code = f"""
int p = {p};
ring R = p, ({','.join(Z)}, t, {avars}, {bvars}), dp;
poly F = {Fstr};
poly Fl = F;
""" + '\n'.join(f'Fl = subst(Fl, {sb});' for sb in subs) + f"""
matrix C = coeffs(Fl, t);
ideal I; int k;
for (k=1;k<=nrows(C);k++) {{ I[k] = C[k,1]; }}
ring S = p, ({avars}, {bvars}), ds;
ideal I = imap(R, I);
ideal G = std(I);
"LOCALDIM"; dim(G);
"LOCALMULT"; mult(G);
quit;
"""
    out = subprocess.run(['Singular', '-q'], input=code, capture_output=True, text=True, timeout=1800).stdout
    toks = out.split()
    d = int(toks[toks.index('LOCALDIM')+1]); m = int(toks[toks.index('LOCALMULT')+1])
    return d, m, out[-300:]
cases = {'Y0=Y1=0 (tdim 4)': ((3,4), 0, roots[0], roots[1]),
         'X0=X1=0 pairing0 (tdim ?)': ((0,1), 0, roots[0], roots[1]),
         'X0=X1=0 pairing1': ((0,1), 1, roots[0], roots[1]),
         'X0=X1=0 pairing2': ((0,1), 2, roots[0], roots[1]),
         'X1=X2=0 (tdim 2)': ((1,2), 0, roots[0], roots[1]),
         'Y0=Y2=0 (tdim 2)': ((3,5), 0, roots[0], roots[1]),
         'X0=X2=0 (tdim 1)': ((0,2), 0, roots[0], roots[1]),
         'X0=Y0=0 pairing0': ((0,3), 0, roots[0], roots[1]),
         'X0=Y0=0 pairing1': ((0,3), 1, roots[0], roots[1]),
         'X0=Y0=0 pairing2': ((0,3), 2, roots[0], roots[1])}
res = {}
for name, (bp, pi, z1, z2) in cases.items():
    pv, qv, (k,l,m,n) = line(bp, pi, z1, z2)
    t0 = time.time()
    try:
        d, mu, tail = singular_local_dim(pv, qv, (l, n))
        res[name] = {'local_dim': d, 'mult': mu, 'line_p': pv, 'line_q': qv, 'sec': round(time.time()-t0,1)}
        print(name, 'local dim', d, 'mult', mu, round(time.time()-t0,1), 's', flush=True)
    except Exception as ex:
        print(name, 'ERR', ex, flush=True)
json.dump(res, open('/tmp/hodge/results/reflection-coupled-residual/construct/lines_local.json','w'), indent=1)
