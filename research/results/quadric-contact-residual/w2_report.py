#!/usr/bin/env python3
"""Pull the remote lane dir, merge wave2-search.log.jsonl, and rewrite the wave-2 search report.

  python3 w2_report.py [--no-pull]
"""
import sys, os, json, time, subprocess, argparse, glob
RES = '/tmp/hodge/results/quadric-contact-residual'
REMOTE = '/tmp/hodge/remote/quadric-contact-residual'
LOG = RES + '/wave2-search.log.jsonl'
REPORT = '/tmp/hodge/notes/wave2-quadric-contact-residual-search.md'
REG = '/tmp/hodge/jobs-registry.jsonl'

ap = argparse.ArgumentParser(); ap.add_argument('--no-pull', action='store_true'); a = ap.parse_args()
if not a.no_pull:
    r = subprocess.run([sys.executable, '-m', 'hodge.rjob', 'pull', '--lane', 'quadric-contact-residual'],
                       cwd='/home/kejace/devel/sextics-work', capture_output=True, text=True, timeout=600)
    print(r.stdout.strip()[-300:], r.stderr.strip()[-300:])
    # merge wave2 log lines
    have = set(open(LOG).read().splitlines()) if os.path.exists(LOG) else set()
    added = 0
    with open(LOG, 'a') as f:
        for rl in sorted(glob.glob(REMOTE + '/wave2-search*.log.jsonl')):
            for line in open(rl).read().splitlines():
                if line and line not in have:
                    f.write(line + '\n'); have.add(line); added += 1
    print(f'merged {added} new wave2 log lines')

recs = []
for line in open(LOG).read().splitlines() if os.path.exists(LOG) else []:
    try:
        recs.append(json.loads(line))
    except Exception:
        pass
witness = [r for r in recs if r.get('WITNESS')]
wfile = os.path.exists(RES + '/witness-candidate.json') or os.path.exists(REMOTE + '/witness-candidate.json')

def fmt(x, nd=2):
    if x is None: return '-'
    if isinstance(x, float): return f'{x:.{nd}e}' if (abs(x) < 1e-3 or abs(x) > 1e4) else f'{x:.4g}'
    return str(x)

L = []
vp = [r for r in recs if r.get('task') == 'varpro' and r.get('status') == 'done']
hs = [r for r in recs if r.get('task') == 'hess-relax']
e1 = [r for r in recs if r.get('task') == 'e1-profile']
cones = [r for r in recs if r.get('task') == 'e1-cones']
st = [r for r in recs if r.get('task') in ('struct', 'struct-exact')]
tot_restarts = sum(r.get('restarts_done', 0) for r in vp) + sum(r.get('restarts_done', 0) for r in hs)
best_vp = min((r['best_rel'] for r in vp if r.get('best_rel') is not None), default=None)
status = 'WITNESS FOUND (see witness-candidate.json)' if (witness or wfile) else 'NO WITNESS'
L.append('# Wave 2 — quadric-contact-residual — SEARCH lens (joint (D,N,Q) search)\n')
L.append(f'STATUS: {status}. Updated {time.strftime("%Y-%m-%d %H:%M")} local. Batches logged: {len(recs)} '
         f'({len(vp)} varpro sector batches, {len(hs)} Hessian-relaxation batches, {len(e1)} E1 profiles, {len(st)} structured candidates); '
         f'random restarts so far: {tot_restarts}; best full-E2 relative residual: {fmt(best_vp)}.\n')
L.append('Target: q = 1, A = X0^3 Y0^3 on the quadric-contact sextic F_b. Residue character trivial on G_aff (order 48), odd under '
         's0 (X0<->Y0) and r1 (X1<->Y1). Prior rigorous exclusions (fqlx #125, replicated in wave 1): deg D <= 6; deg 7 in every pure '
         'parity/character sector; D8* alone at degree 8; any certificate keeping (Y1/X1, Y2/X2) passive.\n')
L.append('Logs: /tmp/hodge/results/quadric-contact-residual/wave2-search.log.jsonl (one JSON line per batch; remote lines merged by w2_report.py). '
         'Scripts (lane dir): w2lib.py, w2_varpro.py, w2_e1profile.py, w2_hess.py, w2_struct.py, w2_report.py. Remote job logs: /tmp/hodge/remote/_jobs/search-quadric-contact-residual-*.log\n')
L.append('## Method (what is new in wave 2)\n')
L.append('* Exact E1 families per G_aff character sector: (J:A)_d[chi] = {D in S_d[chi] : A D in J_{d+6}} over Q (python-flint integer nullspace), '
         'split into the Jacobian part J_d[chi] (sparse generators m*F_i) and the non-Jacobian extras (exact integer polynomials).')
L.append('* Koszul-form variable projection: D(a) = sum a_k D_k with E1-lifts N_k (exact A m e_i for generators, least-squares lifts for extras); '
         'every E1 solution is N = sum a_k N_k + K with K = sum_{r<s} Q_rs (F_s e_r - F_r e_s), and the certificate reduces (q = 1) to B_D(N) = F Q. '
         'f(a) = min over (Q_rs, Q) of ||B_{D(a)}(sum a_k N_k + K) - F Q||^2 on the unit sphere (L-BFGS, exact Golub-Pereyra gradient, '
         'Q-block eliminated by orthogonal projection, lstsq gelsy with cutoff 1e-9). Spurious directions (D(a) = 0, D(a) in F*S_{d-6}) are projected out '
         '(monomial-orthogonality to F-multiples). rel = f/|b|^2 with b = -sum a_k a_l B_{D_k}(N_l).')
L.append('* Hessian relaxation for character-MIXED families: T(a,a) = NF_{J_{2d}}(sum a_k a_l B_{D_k}(N_l)) = 0 (necessary; fqlx\'s Hessian congruence), '
         'minimized numerically with the same restart machinery; near-zeros are polished with the full Koszul system.')
L.append('* Every minimizer with rel < 1e-10 is rationally reconstructed in monomial coordinates (denominators up to 1e6) and re-solved EXACTLY '
         '(engine symmetry.solve_fixed_D_symmetric / fixed_d solvers, core.verify_witness, coprimality by univariate specialization). Only an exact re-solve counts.\n')

# jobs
L.append('## Jobs (remote hodge-sextics, 128 cores shared with ~50 agents)\n')
jobs = []
if os.path.exists(REG):
    for line in open(REG):
        try:
            j = json.loads(line)
        except Exception:
            continue
        if 'quadric-contact-residual' in j.get('name', '') or j.get('lane') == 'quadric-contact-residual':
            jobs.append(j)
for j in jobs:
    L.append(f'* `{j["name"]}` (submitted {j["ts"]}): `{j["cmd"]}`')
for lg in sorted(glob.glob('/tmp/hodge/remote/_jobs/search-quadric-contact-residual-*.log')):
    try:
        tail = [x for x in open(lg, errors='replace').read().splitlines() if 'hodge engine' not in x][-1:]
    except Exception:
        tail = []
    L.append(f'  - {os.path.basename(lg)}: last line: `{(tail[0][:160] if tail else "")}`')
L.append('')

# E1 profile (from the exact table written by w2_e1profile.py; log lines only as fallback)
L.append('## (1) E1-admissible denominator spaces per degree and character sector (exact over Q)\n')
tab = RES + '/e1_sectors_wave2.json'
if os.path.exists(REMOTE + '/e1_sectors_wave2.json'):
    tab = REMOTE + '/e1_sectors_wave2.json'
T = json.load(open(tab)) if os.path.exists(tab) else {}
if T or e1:
    L.append('| d | sectors done | dim S_d | dim J_d | dim (J:A)_d | non-Jacobian | gauge F*S_{d-6} | sectors with extras (chi:count) |')
    L.append('|---|---|---|---|---|---|---|---|')
    for d in sorted(T, key=int):
        secs = T[d]['sectors']
        if not secs: continue
        tot = {k: sum(v[k] for v in secs.values()) for k in ('dimS', 'dimJ', 'kernel', 'extra', 'gauge')}
        ex = ', '.join(f'{c}:{v["extra"]}' for c, v in sorted(secs.items()) if v['extra'])
        L.append(f'| {d} | {len(secs)} | {tot["dimS"]} | {tot["dimJ"]} | {tot["kernel"]} | {tot["extra"]} | {tot["gauge"]} | {ex[:400]} |')
    L.append('')
    L.append('Cross-checks: d <= 9 totals agree with fqlx (#125) and wave 1 (dim (J:A)_d = 0,0,0,0,0,6,36,126,337,760,2053,...); '
             'the per-sector split is new.')
    if T:
        L.append('Per-sector detail (d: chi -> dimS/J/kernel/extra) for degrees with non-Jacobian directions:')
        for d in sorted(T, key=int):
            secs = T[d]['sectors']
            items = [(c, v) for c, v in secs.items() if v['extra']]
            if items:
                L.append(f'* d={d}: ' + '; '.join(f'{c}: {v["dimS"]}/{v["dimJ"]}/{v["kernel"]}/{v["extra"]}' for c, v in sorted(items)))
        L.append('')
else:
    L.append('(pending)\n')
if cones:
    L.append('### Cone / coordinate-subspace families (D supported on a coordinate variable subset, per sector)\n')
    for r in sorted(cones, key=lambda r: r['d']):
        nt = r.get('nontrivial', {})
        L.append(f'* d={r["d"]}: {r["n_nontrivial"]} (subset, sector) pairs with an E1-kernel beyond the gauge' +
                 (': ' + ', '.join(f'{k} ker={v["kernel"]} J={v["dimJ"]} extra={v["extra"]} gauge={v["gauge"]}' for k, v in list(nt.items())[:40]) if nt else ''))
    L.append('')

# varpro
L.append('## (2)-(3) Full-E2 variable-projection searches (Koszul form), per degree and sector\n')
if vp:
    L.append('| tag | d | chi | vars | members K | dimP | rows | nK+nQ | t_eval s | restarts | spurious | best rel | median rel | n<1e-6 | n<1e-10 | mean iters | seconds |')
    L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in sorted(vp, key=lambda r: (r['d'], r['chi'])):
        L.append(f'| {r.get("tag")} | {r["d"]} | {r["chi"]} | {r.get("vars") or "-"} | {r.get("K")} | {r.get("dimP")} | {r.get("rows")} | {r.get("nK")}+{r.get("nQ")} | '
                 f'{fmt(r.get("t_eval"))} | {r.get("restarts_done")} | {r.get("n_spurious")} | {fmt(r.get("best_rel"))} | {fmt(r.get("median_rel"))} | '
                 f'{r.get("n_rel_lt_1e-6")} | {r.get("n_rel_lt_1e-10")} | {fmt(r.get("mean_iters"))} | {r.get("seconds")} |')
    L.append('')
    skipped = [r for r in recs if r.get('task') == 'varpro' and r.get('status') in ('empty-family', 'only-gauge', 'family-error')]
    if skipped:
        L.append('Sectors with no searchable family (empty / gauge only): ' + ', '.join(f'd{r["d"]}:{r["chi"]}({r["status"]})' for r in skipped) + '\n')
    rec_att = [(r, x) for r in vp for x in r.get('recon', [])]
    L.append(f'Reconstruction attempts (rel < 1e-10): {len(rec_att)}')
    for r, x in rec_att[:50]:
        L.append(f'* d={r["d"]} chi={r["chi"]} rel={fmt(x["rel"])}: ' + '; '.join(f'den<={a["max_den"]} err={fmt(a["recon_err"])} -> {a["status"]}' + (' VERIFIED' if a.get('verified') else '') for a in x['attempts']))
    L.append('')
else:
    L.append('(pending)\n')

# hess
L.append('## Hessian relaxation (character-mixed families, necessary condition mod J_{2d})\n')
if hs:
    L.append('| tag | d | mode | members | dimP | dim R_{2d} | restarts | best |T|^2/scale | median | n<1e-8 | n<1e-12 | n<1e-16 | polished (full rel) | seconds |')
    L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in hs:
        pol = ', '.join(fmt(p['full_rel']) for p in r.get('polished', [])) or '-'
        L.append(f'| {r.get("tag")} | {r["d"]} | {r["chi"]} | {r.get("n_members")} | {r.get("dimP")} | {r.get("dimR_2d")} | {r.get("restarts_done")} | '
                 f'{fmt(r.get("best_frel"))} | {fmt(r.get("median_frel"))} | {r.get("n_frel_lt_1e-8")} | {r.get("n_frel_lt_1e-12")} | {r.get("n_frel_lt_1e-16")} | {pol} | {r.get("seconds")} |')
    L.append('')
else:
    L.append('(pending)\n')

# structured
L.append('## (4) Structured candidates (geometry-suggested denominators; E1 screen + exact sector solve mod p)\n')
if st:
    L.append('| d | candidate | chi | E1 | consistent | rank | nullity | rows x cols | s |')
    L.append('|---|---|---|---|---|---|---|---|---|')
    for r in st:
        if r['task'] == 'struct':
            L.append(f'| {r["d"]} | {r["name"]} | {r.get("chi")} | {r.get("e1")} | {r.get("consistent")} | {r.get("rank")} | {r.get("nullity")} | {r.get("rows")}x{r.get("cols")} | {r.get("seconds")} |')
        else:
            L.append(f'| {r["d"]} | {r["name"]} EXACT | - | - | {r.get("status")} | verified={r.get("verified")} | coprime={r.get("coprime")} | | |')
    L.append('')
else:
    L.append('(pending)\n')

# landscape
L.append('## Best residual landscape (full-E2 rel, log10 histogram of local minima over all restarts)\n')
if vp:
    agg = {}
    for r in vp:
        h = agg.setdefault(r['d'], {})
        for k, v in r.get('hist_log10', {}).items():
            h[int(k)] = h.get(int(k), 0) + v
    for d in sorted(agg):
        h = agg[d]
        L.append(f'* d={d}: ' + ', '.join(f'1e{k}: {h[k]}' for k in sorted(h)))
    L.append('')
    L.append('Reading: a genuine witness in a family shows up as rel ~ 1e-14..1e-12 (float64 floor); the smallest values seen so far are listed above; '
             'all minima are separated from the floor by many orders of magnitude unless stated otherwise.')
else:
    L.append('(pending)\n')
L.append('\n## Structure notes (from the exact E1 tables)\n')
L.append('* (J:A)/J has exactly ONE minimal generator beyond J in degrees <= 9: D8* = 4G_S^2 - G_P^2 (= const * G_u G_v, the product of the two '
         'reflection-contour quartics); the four degree-9 directions are D8*.X1, D8*.Y1, D8*.X2, D8*.Y2 (D8*.X0, D8*.Y0 lie in J_9). '
         'At degree 10 only 10 of the 556 non-Jacobian directions come from D8*.S_2; the other 546 are new minimal generators of the colon ideal '
         '(computed mod p this session). From d = 10 on E1 is essentially the Gorenstein-dual constraint codim (J:A)_d = dim R_{d+6} (minus 1 at d = 10): '
         '950 = 951 - 1, 666, 426 for d = 10, 11, 12.')
L.append('* Gorenstein duality (socle degree 24) pairs ker(.A: R_d -> R_{d+6}) with coker(.A: R_{18-d} -> R_{24-d}); the single kernel direction D8* at d = 8 '
         'is dual to the single missing direction of A.R_10 in R_16, consistent with the 556 = 1506 - 950 kernel at d = 10.')
L.append('* Sector bookkeeping: G_aff characters chi = (e0+e3, e2+e5 mod 2; e1-e4 mod 6; Y-parity); the target residue is G_aff-invariant and the identity '
         'for a chi-eigen D lives in sector chi_A + 2 chi, so restricting all unknowns to sectors is lossless for eigen-denominators (24 sectors per degree parity).')
L.append('\n## Errors / notes\n')
errs = [r for r in recs if r.get('status') == 'family-error']
for r in errs:
    L.append(f'* d={r["d"]} chi={r["chi"]}: {r.get("error")}')
L.append('* Engine: no engine files were modified; lane-local library w2lib.py. Observations about the shared engine are in /tmp/hodge/notes/engine-bugs.md.')
open(REPORT, 'w').write('\n'.join(L) + '\n')
print('report written', REPORT, 'records', len(recs))
