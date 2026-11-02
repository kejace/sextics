"""Regenerate /tmp/hodge/notes/wave2-klein-coupled-character-search.md from the batch logs."""
import json, os, time, collections
LANE_DIR = '/tmp/hodge/results/klein-coupled-character'
OUT = '/tmp/hodge/notes/wave2-klein-coupled-character-search.md'
def read(path):
    rows = []
    if os.path.exists(path):
        for line in open(path):
            try: rows.append(json.loads(line))
            except Exception: pass
    return rows
batches = read(os.path.join(LANE_DIR, 'wave2-search.log.jsonl'))
dims = read(os.path.join(LANE_DIR, 'w2_dims.jsonl'))
search = [b for b in batches if b.get('lens') == 'search' and 'restarts' in b]
witness = any(b.get('WITNESS') for b in batches) or os.path.exists(os.path.join(LANE_DIR, 'witness-candidate.json'))
best = min((b['best_rel'] for b in search if b.get('best_rel') is not None), default=None)
nrest = sum(b.get('restarts', 0) for b in search)
confirms = [(b, p) for b in search for p in b.get('polished', []) if p.get('confirm')]
L = []
L.append('# Wave 2 — klein-coupled-character — SEARCH lens (joint denominator/numerator search)\n')
status = 'WITNESS FOUND' if witness else 'NO WITNESS'
L.append(f'STATUS: {status}. {len(search)} search batches, {nrest} restarts total, best relative residual {best if best is None else format(best, ".3e")}; '
         f'{len(confirms)} exact-confirmation attempts. Updated {time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())}. '
         f'Lane dir /tmp/hodge/results/klein-coupled-character (w2fast.py engine, w2_search.py driver, w2_dims.py, w2_validate.py; '
         f'batch log wave2-search.log.jsonl; per-restart streams wave2-<tag>-d*p*c*-w*.jsonl).\n')
L.append('## Method\n')
L.append('Target q = 1, A = X0^3 Y0^3 (Y-odd, sigma-character 3; weights W = (1,2,4,0,0,0) mod 7, F character 0). '
         'For a sigma-eigen denominator D of degree d, Y-parity parD, character chi the certificate system is restricted losslessly '
         '(averaging) to N_i (deg d+1, parity 1+parD+pi_i, char 3+chi+W_i), Koszul corrections p_ij (deg d-4, parity 1+parD+pi_i+pi_j, '
         'char 3+chi+W_i+W_j), Q (deg 2d-6, Y-odd, char 3+2chi). Correction-map form (ch. 5.2): with an E1 lift N0 (N0(F) = A D) the '
         'condition is B_D(N0) + sum p_ij L_ij = F Q. Joint search = variable projection over a linear family D(a) (E1-admissible '
         'sector space, lifts linear in a): f(a) = min_{p,Q} |T(a) p - F Q - b(a)|^2 on the unit sphere, F-multiples projected out, '
         'L-BFGS random restarts; rel = f/|b|^2. Inner solve: projected normal equations with eigen-cutoff 1e-7 (checked against LAPACK gelsd to 1e-16). '
         'Engine validated on accepted q=1 witnesses (f = 0 at the witness; Koszul-perturbed lifts absorbed to 1e-23 relative), gradient checked by finite differences. '
         'Families: full = (J:A)_{d,sector} (float null space of the sector E1 system = J_d + non-Jacobian kernel); w = J_d with lifts A w; '
         'hess = Hessian relaxation (N = A w exactly, only Q free). A numerical residual is NOT evidence; only exact re-solves (kcc2 Koszul test mod p + multi-prime lift + core.verify_witness) count.\n')
if dims:
    L.append('## E1-admissible denominator spaces per sector (exact mod p = 2147483647; K_float = float null-space dimension cross-check)\n')
    L.append('| d | parD | chi | dim S_d sector | dim J_d | dim (J:A)_d | non-Jacobian kernel | F-multiples | nN | nP | nQ | rows |')
    L.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for r in dims:
        kf = r.get('K_float'); flag = '' if kf is None or kf == r['kernel_total'] else f' (float {kf}!)'
        L.append(f"| {r['d']} | {r['parD']} | {r['chi']} | {r['nD']} | {r['J_d']} | {r['kernel_total']}{flag} | {r['nontrivial']} | {r.get('F_mult','')} | {r['nN']} | {r['nP']} | {r['nQ']} | {r['rows']} |")
    tot = collections.defaultdict(lambda: [0, 0, 0])
    for r in dims:
        t = tot[(r['d'], r['parD'])]; t[0] += r['J_d']; t[1] += r['kernel_total']; t[2] += r['nontrivial']
    L.append('\nTotals per (d, parity): ' + '; '.join(f"d={k[0]} par={k[1]}: J={v[0]}, (J:A)={v[1]}, kernel mod J={v[2]}" for k, v in sorted(tot.items())) + '\n')
L.append('## Search batches\n')
L.append('| family | d | parD | chi | K (dimP) | nP/nQ/rows | restarts (spurious) | best rel | residual histogram | distinct minima (rel, rounded) | rank T at best (nullity) | polish -> rel | s/restart |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for b in search:
    h = ' '.join(f"{k}:{v}" for k, v in sorted(b.get('hist_rel', {}).items(), key=lambda kv: float(kv[0].strip('<>=')) if kv[0][0] == '<' else 99))
    pol = ', '.join(f"{p['rel_after']:.2e}" for p in b.get('polished', []))
    rk = b.get('rankT_at_best'); nul = (b['nP'] - rk) if (rk is not None and b.get('nP')) else ''
    L.append(f"| {b['family']} | {b['d']} | {b['parD']} | {b['chi']} | {b['K']} ({b['dimP']}) | {b['nP']}/{b['nQ']}/{b['rows']} | {b['restarts']} ({b['spurious']}) | "
             f"{b['best_rel']:.3e} | {h} | {', '.join(f'{x:.4g}' for x in b.get('distinct_minima', [])[:6])} | {rk} ({nul}) | {pol} | {b.get('mean_restart_seconds', 0):.0f} |")
empties = [b for b in batches if b.get('lens') == 'search' and 'restarts' not in b]
if empties:
    L.append('\nEmpty / trivial sectors: ' + '; '.join(f"d={b['d']} par={b['parD']} chi={b['chi']} ({b.get('note')})" for b in empties) + '\n')
if confirms:
    L.append('\n## Exact confirmation attempts (rational reconstruction + exact Koszul test)\n')
    for b, p in confirms:
        for c in p['confirm']:
            L.append(f"- d={b['d']} par={b['parD']} chi={b['chi']} family={b['family']}: rel {p['rel_after']:.2e} -> {c['status']} (rank {c.get('rank')}, nullity {c.get('nullity')}); D = {c['D'][:200]}")
other = [b for b in batches if b.get('lens') not in ('search',)]
if other:
    L.append('\n## Other batches (structured families / sweeps / notes)\n')
    for b in other:
        L.append('- ' + json.dumps({k: v for k, v in b.items() if k not in ('time', 'engine', 'instance')}, default=str)[:600])
L.append('\n## Interpretation / caveats\n')
L.append('- A local minimum with rel > 1e-6 says nothing rigorous; a rel ~ 1e-1 .. 2e-1 plateau reached from many restarts suggests no solution in that sector at this degree but is not an exclusion.')
L.append('- Y-even and Y-odd D are both covered; non-eigen D (mixing characters) are not covered by sector-restricted searches.')
L.append('- Exclusions of degree <= 8 outside J(F) and the Jhat / polar supports are rigorous prior results (dossier, attack1, verify lens); this lens adds numerical evidence only unless a witness is confirmed exactly.')
open(OUT, 'w').write('\n'.join(L) + '\n')
print('report written', OUT, len(search), 'batches')
