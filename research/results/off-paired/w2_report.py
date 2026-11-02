"""Rewrite /tmp/hodge/notes/wave2-off-paired-search.md from the lane logs (STATUS first)."""
import json, os, time, glob
LOG = '/tmp/hodge/results/off-paired/wave2-search.log.jsonl'
DIMS = '/tmp/hodge/results/off-paired/wave2-dims.jsonl'
CAND = '/tmp/hodge/results/off-paired/wave2-candidates.jsonl'
OUT = '/tmp/hodge/notes/wave2-off-paired-search.md'


def read(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except Exception:
                    pass
    return out


recs = read(LOG); dims = read(DIMS); cands = read(CAND)
wit = [r for r in recs if r.get('step') == 'WITNESS']
vp = [r for r in recs if r.get('step') == 'varpro' and r.get('restarts', 0) > 0]
sw = [r for r in recs if r.get('step') == 'sweep']
other = [r for r in recs if r.get('step') not in ('varpro', 'WITNESS', 'dims', 'sweep')]
tot_restarts = sum(r['restarts'] for r in vp)
best = min((r['best_rel'] for r in vp), default=None)
L = []
if wit:
    L.append(f"# WITNESS FOUND — off-paired — see /tmp/hodge/results/off-paired/witness-candidate.json ({wit[-1].get('D','')[:80]})")
    L.append('')
L.append('# Wave 2 — off-paired — SEARCH lens')
L.append('')
status = 'WITNESS FOUND' if wit else 'No witness.'
L.append(f"STATUS ({time.strftime('%Y-%m-%d %H:%M')} UTC): {status} {len(vp)} variable-projection batches, "
         f"{tot_restarts} restarts total, best relative residual {best:.2e}" if best is not None else
         f"STATUS ({time.strftime('%Y-%m-%d %H:%M')} UTC): {status} no search batches logged yet.")
L.append(f"Candidates with rel < 1e-7 logged: {len(cands)} (wave2-candidates.jsonl); exact confirmations: {len(wit)}.")
L.append('')
L.append('Lane dir /tmp/hodge/results/off-paired: library w2lib.py (sector-restricted reduced variable projection: D in an '
         'E1-admissible sector basis with Euler gauge removed, N = E1-lift(a) + Koszul unknowns, equation B_D(N) = F Q), '
         'drivers w2_dims.py / w2_varpro.py / w2_sweep_*.py, report generator w2_report.py. Batch log wave2-search.log.jsonl; '
         'per-restart records wave2-<job>.jsonl; dims wave2-dims.jsonl. Engine v2. rjob cannot ssh from the VM to itself, so jobs '
         'are launched with nohup setsid (logs /tmp/hodge/jobs/search-off-paired-*.log).')
L.append('')
L.append('## Method')
L.append('* Characters (g,h) of the diagonal group Z/6 x Z/6, g=(0,2,0,1,3,3), h=(0,0,1,0,0,5); chi(A)=(3,0), chi(F)=(0,0). '
         'For D of character chi: chi(N_i) = chi_A + chi + w_i, chi(Q) = chi_A + 2 chi; rows = degree-2d monomials of character chi_A + 2 chi.')
L.append('* "jac" family = J(F)_{d,chi} modulo F*R_{d-6} (basis m F_i, exact lift N = A m e_i / q); "ker" = complement of J in (J:A)_{d,chi} '
         '(non-Jacobian E1 kernel, float least-squares lifts); "full" = both. "hess" = Hessian relaxation (no Koszul unknowns, N = A w/q).')
L.append('* Optimizer: hodge.varpro L-BFGS on the unit sphere with the Golub-Pereyra gradient; rel = |r|^2/|b|^2 at the minimizer; '
         'rel < 1e-5 polished by trust-region least squares; rel < 1e-10 -> rational reconstruction (denominators to 1e5) + exact '
         'sector-restricted re-solve (symmetry.solve_fixed_D_symmetric). Only an exact re-solve counts.')
L.append('')
# dims table
if dims:
    L.append('## 1. E1-admissible denominator spaces (exact over Q; dims after/before gauge)')
    L.append('')
    bydeg = {}
    for s in dims:
        bydeg.setdefault(s['d'], []).append(s)
    L.append('| d | sectors (E1 > gauge) | sum dim (J:A)_d | sum gauge | sum J mod gauge | sum non-Jacobian kernel | largest sector (chi: E1/gauge/jac/ker) |')
    L.append('|---|---|---|---|---|---|---|')
    for d in sorted(bydeg):
        ss = bydeg[d]
        nz = [s for s in ss if s['dimE1'] - s['gauge'] > 0]
        big = max(ss, key=lambda s: s['dimE1'])
        L.append(f"| {d} | {len(nz)} | {sum(s['dimE1'] for s in ss)} | {sum(s['gauge'] for s in ss)} | "
                 f"{sum(s['jac_mod_gauge'] for s in ss)} | {sum(s['ker'] for s in ss)} | "
                 f"{tuple(big['chi'])}: {big['dimE1']}/{big['gauge']}/{big['jac_mod_gauge']}/{big['ker']} |")
    L.append('')
    L.append('Per-sector detail (d, chi, R_d^chi, E1, gauge, jac mod gauge, ker):')
    L.append('')
    for d in sorted(bydeg):
        row = ', '.join(f"{tuple(s['chi'])}:{s['dimR']}/{s['dimE1']}/{s['gauge']}/{s['jac_mod_gauge']}/{s['ker']}"
                        for s in sorted(bydeg[d], key=lambda s: tuple(s['chi'])) if s['dimE1'] - s['gauge'] > 0)
        L.append(f"* d={d}: {row}")
    L.append('')
if vp:
    L.append('## 2/3. Variable-projection batches (one line per job)')
    L.append('')
    L.append('| job | d | chi | family | hess | K (D-dim) | S (Koszul) | Q | rows | rank M (random a) | restarts | best rel | median rel | log10 q10/q90 | #<1e-4 | #<1e-7 | wall s |')
    L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in sorted(vp, key=lambda r: (r['d'], tuple(r['chi']), r['kind'], r['hess'], r['ts'])):
        fi = r.get('family', {})
        L.append(f"| {r['name']} | {r['d']} | {tuple(r['chi'])} | {r['kind']} | {r['hess']} | {fi.get('K','')} | {fi.get('nS','')} | {fi.get('nQ','')} | "
                 f"{fi.get('rows','')} | {fi.get('rank_M_random','')} | {r['restarts']} | {r['best_rel']:.2e} | {r['median_rel']:.2e} | "
                 f"{r.get('q10_log10',0):.1f}/{r.get('q90_log10',0):.1f} | {r.get('n_below_1e-4',0)} | {r.get('n_below_1e-7',0)} | {r['seconds_wall']} |")
    L.append('')
if sw:
    L.append('## 4. Structured-family sweeps')
    L.append('')
    L.append('| batch | family | degrees | candidates | E1 pass | consistent | note |')
    L.append('|---|---|---|---|---|---|---|')
    for r in sw:
        L.append(f"| {r.get('name','')} | {r.get('family','')} | {r.get('degrees','')} | {r.get('n',0)} | {r.get('e1_pass',0)} | {r.get('consistent',0)} | {r.get('note','')} |")
    L.append('')
if cands:
    L.append('## 5. Candidates (rel < 1e-7) and exact re-solves')
    L.append('')
    for c in cands[-30:]:
        L.append(f"* {c.get('name')} d={c.get('d')} chi={c.get('chi')} rel={c.get('rel'):.2e} polished={c.get('polished_rel','-')} D~ {c.get('D_float','')[:140]}")
    L.append('')
if other:
    L.append('## Other log entries')
    for r in other[-20:]:
        L.append('* ' + json.dumps({k: v for k, v in r.items() if k not in ('family',)})[:400])
    L.append('')
L.append('## Notes / interpretation')
L.append('* A stalled residual is not evidence of nonexistence; only exact exclusions (msolve/Singular, degrees 6-7 all sectors, '
         'degree 8 outside the six (a,0) sectors — attack1 and deanodotdotdot) are theorems-mod-p.')
L.append('* Residual floors: jac-type families on off-paired typically stall at rel ~1e-2..1e-4 (see table); a genuine zero would show '
         'rel < 1e-12 after polishing (Fermat control: D = f6 found with rel ~1e-30 and confirmed exactly).')
tmp = OUT + '.tmp'
with open(tmp, 'w') as f:
    f.write('\n'.join(L) + '\n')
os.replace(tmp, OUT)
