"""R2: exact facts about D8 and exact fixed-D exclusions.
(a) the packet's D8 (mirror e1_d8_exact.json): sector, D8 not in J_8, A*D8 in J_14, sum_i P_i d_i(7F) = A*D8 exactly;
(b) own exact D8 from an fmpz nullspace, congruent to the packet's mod J_8;
(c) own single-identity (formal contract) G-sector system for D8: exact fmpz ranks;
(d) replay of the parallel session's exact Farkas duals (receipts/farkas_*.json) against ALL columns, no parity restriction."""
import sys, json, time
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/results/cyclic-triple-residual/wave2-referee')
import ref_lib as L
from flint import fmpz_mat
MIR = '/tmp/hodge/results/cyclic-triple-residual/wave2-referee/mirror'   # copies of the remote-only receipts (receipts/, e1_d8_exact.json)
OUT = '/tmp/hodge/results/cyclic-triple-residual/wave2-referee'
LOG = OUT + '/referee.log.jsonl'
def log(rec):
    rec['ts'] = time.strftime('%Y-%m-%dT%H:%M:%S'); rec['script'] = 'ref_d8.py'
    with open(LOG, 'a') as f: f.write(json.dumps(rec, default=str) + '\n')
L.check_setup()
t0 = time.time()
dat = json.load(open(MIR + '/e1_d8_exact.json'))
D8 = {tuple(map(int, k.split(','))): v for k, v in dat['D8z'].items()}
Pz = [{tuple(map(int, k.split(','))): v for k, v in Pi.items()} for Pi in dat['Pz']]
cD = {L.chi(e) for e in D8}
print('D8 terms', len(D8), 'max|c|', max(abs(v) for v in D8.values()), 'characters', cD, flush=True)
inJ8 = L.in_J_exact(D8)
AD8 = L.pmul(L.AZ, D8)
inJ14 = L.in_J_exact(AD8)
lhs = {}
for i in range(6): lhs = L.padd(lhs, L.pmul(Pz[i], L.DF[i]))
pid = (lhs == AD8)
print(f'D8 in J_8: {inJ8}; A*D8 in J_14: {inJ14}; sum P_i d_i(7F) == A*D8: {pid}  ({time.time()-t0:.1f}s)', flush=True)
log({'check': 'D8_exact_facts', 'terms': len(D8), 'sector': str(sorted(cD)), 'D8_in_J8': inJ8[0], 'ranks_J8': inJ8[1:], 'AD8_in_J14': inJ14[0],
     'ranks_J14': inJ14[1:], 'P_identity_exact': pid, 'seconds': round(time.time() - t0, 1)})

# (b) own exact kernel at d = 8, sector (0,4): nullspace of [J_14 gens | A*m] over Z
t1 = time.time()
c = (0, 4); c2 = L.cadd(c, L.CHI_A)
S8 = L.mons_sector(8, c); tgt = L.mons_sector(14, c2); idx = {m: i for i, m in enumerate(tgt)}
gens = L.jac_gens(14, c2); cols = [L.pmul(L.AZ, L.mono(m)) for m in S8]
allc = gens + cols
flat = L.to_rows(allc, idx)
M = fmpz_mat(len(allc), len(tgt), flat).transpose()          # columns = generators and A*m
X, nullity = M.nullspace()
# kernel D-parts: last len(S8) coordinates of each null vector
Dparts = []
for k in range(nullity):
    v = [int(X[len(gens) + j, k]) for j in range(len(S8))]
    if any(v): Dparts.append({S8[j]: v[j] for j in range(len(S8)) if v[j]})
idx8 = {m: i for i, m in enumerate(S8)}
J8 = L.jac_gens(8, c)
rJ8 = L.rank_rows(J8, idx8); rK = L.rank_rows(J8 + Dparts, idx8); rKD = L.rank_rows(J8 + Dparts + [D8], idx8)
print(f'own kernel: nullity {nullity}, D-parts span rank {rK} vs J_8 rank {rJ8} => Ann_8 = {rK - rJ8}; adding packet D8: rank {rKD} ({time.time()-t1:.1f}s)', flush=True)
log({'check': 'D8_own_exact_kernel', 'sector': str(c), 'nullity_fmpz': nullity, 'rank_J8': rJ8, 'rank_J8_plus_kernel': rK, 'ann8_exact': rK - rJ8,
     'rank_with_packet_D8': rKD, 'packet_D8_in_span': rKD == rK, 'seconds': round(time.time() - t1, 1)})

# (c) own fixed-D single-identity sector system for D8, exact ranks
t2 = time.time()
cols, meta, rhs = L.fixedD_columns(D8, sector_only=True)
nN = sum(1 for m in meta if m[0] == 'N'); nQ = len(meta) - nN
res = L.exact_ranks(cols, rhs)
print(f'own fixed-D system for D8 (single identity, G-sector): N cols {nN}, Q cols {nQ}, {res} ({time.time()-t2:.1f}s)', flush=True)
log({'check': 'D8_fixedD_exact', 'identity': 'single cleared identity E = F^2 Q (formal contract), G-sector reduced', 'N_cols': nN, 'Q_cols': nQ, **res,
     'seconds': round(time.time() - t2, 1)})

# (d) Farkas replays
def parse_D(s):
    # tiny parser for the receipt's D strings: terms like  -54*X0**6*Y0*Y1/5  (sums of monomials with rational coeffs)
    import re
    s = s.replace(' ', '')
    terms = re.findall(r'[+-]?[^+-]+', s)
    P = {}
    names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
    for t in terms:
        sign = -1 if t.startswith('-') else 1
        t = t.lstrip('+-')
        den = 1
        if '/' in t: t, den = t.split('/'); den = int(den)
        parts = t.split('*')
        coef = Fraction(1)
        e = [0] * 6
        i = 0
        while i < len(parts):
            q = parts[i]
            if q in names:
                k = 1
                if i + 1 < len(parts) and parts[i + 1] == '':   # '**' produced an empty split element
                    k = int(parts[i + 2]); i += 2
                e[names.index(q)] += k
            else:
                coef *= Fraction(int(q))
            i += 1
        coef = coef * sign / den
        e = tuple(e)
        P[e] = P.get(e, 0) + coef
    from math import lcm
    Lc = 1
    for v in P.values(): Lc = lcm(Lc, Fraction(v).denominator)
    return {e: int(v * Lc) for e, v in P.items() if v}

for name in ['farkas_D8', 'farkas_E_X_F', 'farkas_fqlx_cone1', 'farkas_fqlx_cone2']:
    t3 = time.time()
    rec = json.load(open(f'{MIR}/receipts/{name}.json'))
    if not rec.get('farkas'):
        print(name, 'no Farkas vector in receipt'); continue
    Dz = parse_D(rec['D'])
    y = {tuple(int(x) for x in k.strip('()').split(',')): Fraction(v) for k, v in rec['farkas']['vector'].items()}
    cols, meta, rhs = L.fixedD_columns(Dz, sector_only=False)      # ALL columns, no parity / sector restriction
    bad = 0
    for col in cols:
        s = sum(v * y[e] for e, v in col.items() if e in y)
        if s != 0: bad += 1
    val = sum(v * y[e] for e, v in rhs.items() if e in y)
    ok = (bad == 0 and val != 0)
    print(f'{name}: D terms {len(Dz)} deg {L.pdeg(Dz)}; y nonzeros {len(y)}; all {len(cols)} columns orthogonal: bad={bad}; y.rhs={val}; '
          f'Farkas certificate valid: {ok} ({time.time()-t3:.1f}s)', flush=True)
    log({'check': f'farkas_replay_{name}', 'D_terms': len(Dz), 'deg': L.pdeg(Dz), 'y_nonzeros': len(y), 'columns_checked': len(cols), 'bad_columns': bad,
         'y_dot_rhs': str(val), 'valid': ok, 'seconds': round(time.time() - t3, 1)})
print('done', flush=True)
