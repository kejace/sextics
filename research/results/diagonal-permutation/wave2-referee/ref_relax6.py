"""Independent construction of the degree-6 Hessian relaxation  A * H(w,w) in J_12  for linear vector fields w,
sectors e (D Y-even), o (D Y-odd), m (mixed), at a fresh prime; Euler gauge: coefficient of X0 in w_{X0} set to 0.
Writes ref_relax6<s>_p<p>.ms, runs msolve -g 2, checks pure powers of every variable among the leading terms
(grevlex), counts standard monomials (colength).  Receipt ref_relax6_p<p>.json."""
import sys, json, time, subprocess, re, itertools; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation/wave2-referee')
import ref_common as R
from flint import nmod_mat
p = int(sys.argv[1]) if len(sys.argv) > 1 else 999999937
sectors = sys.argv[2].split(',') if len(sys.argv) > 2 else ['e', 'o', 'm']
F, A, q = R.load_instance()
plus, minus = R.koszul_series(30)
parts = [R.pderiv(F, i) for i in range(6)]
hess = {(i, j): {e: R.modp(c, p) for e, c in R.pderiv(parts[i], j).items()} for i in range(6) for j in range(6)}
Ap = {e: R.modp(c, p) for e, c in A.items()}
partsp = [{e: R.modp(c, p) for e, c in P.items()} for P in parts]
lin = R.monomials(1)
ypar = lambda e: (e[3] + e[4] + e[5]) % 2

def jac_columns(n, par):
    cols = []
    for i in range(6):
        mpar = par if i < 3 else 1 - par
        for m in R.monomials(n - 5, mpar):
            cols.append(R.mul({m: 1}, partsp[i]))
    return cols

def left_nullspace(rows_mons, cols):
    """basis of functionals on S_n^par vanishing on span(cols): returns (L as nmod_mat with dim R rows, idx)."""
    idx = {m: k for k, m in enumerate(rows_mons)}
    Mt = nmod_mat(len(cols), len(rows_mons), p)        # transposed generator matrix
    for j, P in enumerate(cols):
        for e, c in P.items():
            Mt[j, idx[e]] = c % p
    X, nullity = Mt.nullspace()
    # columns 0..nullity-1 of X span the nullspace of Mt, i.e. functionals killing all generators
    L = nmod_mat(nullity, len(rows_mons), p)
    for a in range(nullity):
        for b in range(len(rows_mons)):
            v = int(X[b, a])
            if v: L[a, b] = v
    return L, idx, nullity

t0 = time.time()
NF = {}
for par in (0, 1):
    rows = R.monomials(12, par)
    cols = jac_columns(12, par)
    L, idx, nullity = left_nullspace(rows, cols)
    good = (nullity == (plus if par == 0 else minus)[12])
    NF[par] = (L, idx, nullity, good)
    print(f'degree-12 parity {par}: dim S = {len(rows)}, generators {len(cols)}, dim R_12 = {nullity} (Koszul {(plus if par==0 else minus)[12]}) good_reduction={good}', flush=True)
receipt = {'p': p, 'good_reduction_12': {par: NF[par][3] for par in (0, 1)}, 'sectors': {}}

# unknowns: (i, m) with w_i = sum_m c_{i,m} m ; D-parity of the term c_{i,m} m F_i is ypar(m) + (0 if i<3 else 1)
unknowns_all = [(i, m) for i in range(6) for m in lin]
def dpar(u):
    i, m = u; return (ypar(m) + (0 if i < 3 else 1)) % 2
for s in sectors:
    eps = {'e': [0], 'o': [1], 'm': [0, 1]}[s]
    unk = [u for u in unknowns_all if dpar(u) in eps and not (u[0] == 0 and u[1] == (1, 0, 0, 0, 0, 0))]  # Euler gauge
    names = [f'c{i}_{"".join(map(str, m))}' for i, m in unk]
    n = len(unk)
    pairs = [(a, b) for a in range(n) for b in range(a, n)]
    # Phi_ab = A * m_a m_b F_{i_a i_b} * (2 if a != b else 1)  -- a polynomial of degree 12 with definite parity
    cols_by_par = {0: [], 1: []}; pair_par = []
    for (a, b) in pairs:
        (ia, ma), (ib, mb) = unk[a], unk[b]
        P = R.mul(R.mul({ma: 1}, {mb: 1}), hess[(ia, ib)])
        if a != b: P = {e: (2 * c) % p for e, c in P.items()}
        P = R.mul(Ap, P)
        P = {e: c % p for e, c in P.items() if c % p}
        par = (ypar(ma) + ypar(mb) + (ia >= 3) + (ib >= 3) + 1) % 2  # parity of A * m_a m_b F_{ia ib}
        pair_par.append(par); cols_by_par[par].append((len(pair_par) - 1, P))
    eqs = []
    for par in (0, 1):
        L, idx, nullity, good = NF[par]
        if not cols_by_par[par]: continue
        Phi = nmod_mat(len(idx), len(pairs), p)
        for (k, P) in cols_by_par[par]:
            for e, c in P.items():
                Phi[idx[e], k] = c
        E = L * Phi                       # nullity x pairs : coordinates of A*Phi_ab in R_12^par
        Er, rk = E.rref()
        for r in range(rk):
            eqs.append([int(Er[r, k]) for k in range(len(pairs))])
    Eall = nmod_mat(eqs, p) if eqs else None
    Eall_r, rank_all = Eall.rref()
    rows_ind = [[int(Eall_r[r, k]) for k in range(len(pairs))] for r in range(rank_all)]
    ms = f'ref_relax6{s}_p{p}.ms'
    with open(ms, 'w') as f:
        f.write(','.join(names) + '\n' + str(p) + '\n')
        polys = []
        for row in rows_ind:
            terms = []
            for k, (a, b) in enumerate(pairs):
                if row[k]:
                    mon = f'{names[a]}^2' if a == b else f'{names[a]}*{names[b]}'
                    terms.append(f'{row[k]}*{mon}')
            polys.append('+'.join(terms))
        f.write(',\n'.join(polys) + '\n')
    gb = f'ref_relax6{s}_p{p}.gb'
    t1 = time.time()
    pr = subprocess.run(['msolve', '-f', ms, '-o', gb, '-t', '2', '-g', '2', '-v', '1'], capture_output=True, text=True, timeout=3600)
    txt = open(gb).read()
    body = '\n'.join(l for l in txt.splitlines() if not l.startswith('#')).strip().lstrip('[').rstrip(':').rstrip(']')
    gpolys = [x.strip() for x in body.replace('\n', '').split(',') if x.strip()]
    nidx = {nm: i for i, nm in enumerate(names)}
    def lead(poly):
        best = None
        for term in re.split(r'(?=[+-])', poly.replace(' ', '')):
            if not term or term in '+-': continue
            expo = [0] * n
            for fac in term.lstrip('+-').split('*'):
                if not fac or re.fullmatch(r'\d+', fac): continue
                if '^' in fac:
                    v, e = fac.split('^'); expo[nidx[v]] += int(e)
                else: expo[nidx[fac]] += 1
            key = (sum(expo), tuple(-x for x in reversed(expo)))
            if best is None or key > best[0]: best = (key, tuple(expo))
        return best[1]
    leads = [lead(g) for g in gpolys]
    pure = {}
    for l in leads:
        nz = [i for i, e in enumerate(l) if e]
        if len(nz) == 1: pure[nz[0]] = min(pure.get(nz[0], 10**9), l[nz[0]])
    zero_dim = all(i in pure for i in range(n))
    colength = None
    if zero_dim:
        bounds = [pure[i] for i in range(n)]
        cnt = 0
        for e in itertools.product(*[range(b) for b in bounds]):
            if not any(all(e[i] >= l[i] for i in range(n)) for l in leads): cnt += 1
        colength = cnt
    rec = {'unknowns': n, 'pair_monomials': len(pairs), 'independent_equations': rank_all, 'gb_size': len(gpolys),
           'gb_max_lead_degree': max(sum(l) for l in leads), 'pure_powers': len(pure), 'zero_dimensional_V_is_origin': zero_dim,
           'colength_of_leading_ideal': colength, 'msolve_seconds': round(time.time() - t1, 2), 'msolve_status': pr.returncode,
           'dim_quotient_reported_by_msolve': (re.findall(r'Dimension of quotient:\s*(\d+)', pr.stdout + pr.stderr) or [None])[-1]}
    receipt['sectors'][s] = rec
    print(s, json.dumps(rec), flush=True)
receipt['seconds'] = round(time.time() - t0, 1)
json.dump(receipt, open(f'ref_relax6_p{p}.json', 'w'), indent=1)
print('DONE')
