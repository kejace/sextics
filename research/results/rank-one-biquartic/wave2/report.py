"""Rebuild /tmp/hodge/notes/wave2-rank-one-biquartic-search.md from the merged logs (STATUS first)."""
import json, os, glob, time, statistics
LANE = '/tmp/hodge/results/rank-one-biquartic'
MIR = '/tmp/hodge/remote/rank-one-biquartic/wave2/out'
LOG = f'{LANE}/wave2-search.log.jsonl'
REP = '/tmp/hodge/notes/wave2-rank-one-biquartic-search.md'
recs = [json.loads(l) for l in open(LOG) if l.strip()]
def _is_lane_witness(r):
    tag = str(r.get('tag', ''))
    if tag.startswith(('planted', 'smoke')) or r.get('instance') not in (None, 'rank-one-biquartic'):
        return False
    return bool(r.get('witness') or r.get('WITNESS') or r.get('status') == 'WITNESS')
wit = any(_is_lane_witness(r) for r in recs)
wit = wit or os.path.exists(f'{LANE}/witness-candidate.json') or os.path.exists('/tmp/hodge/remote/rank-one-biquartic/witness-candidate.json')
out = []
out.append('# Wave 2 search report: `rank-one-biquartic` (q = 1, A = X0^3 Y0^3)\n')
out.append(('**STATUS: WITNESS FOUND -- see section 6 and witness-candidate.json**' if wit else
            '**STATUS: NO WITNESS FOUND** (extensive wave-2 search; details below; search log `wave2-search.log.jsonl`).') + '\n')
out.append(f'Last rebuilt {time.strftime("%Y-%m-%d %H:%M:%S")}.  Work dir `/tmp/hodge/results/rank-one-biquartic/wave2/` '
           '(scripts; `out/` job outputs on the remote, mirrored to `/tmp/hodge/remote/rank-one-biquartic/wave2/out/`).\n')
# static method section
out.append(open(os.path.join(os.path.dirname(__file__), 'report_method.md')).read())
# 1. E1 spaces
e1 = sorted([r for r in recs if r.get('batch') == 'e1spaces'], key=lambda r: (r['d'], r['parity']))
out.append('\n## 1. E1-admissible denominator spaces (exact mod p = 2^31-1)\n')
if e1:
    out.append('| d | parity | dim S_d | dim J_d | ring dim (S/J)_d | dim {D: AD in J} | kernel mod J | matrix | time s |\n|---|---|---|---|---|---|---|---|---|\n')
    seen = set()
    for r in e1:
        k = (r['d'], r['parity'])
        if k in seen: continue
        seen.add(k)
        out.append(f"| {r['d']} | {'even' if r['parity']==0 else 'odd'} | {r['nD']} | {r['J_dim']} | {r.get('ring_dim','')} | {r['e1_dim']} | **{r['kernel_mod_J']}** | {r['rows']}x{r['cols']} | {r['t_build']+r['t_null']:.0f} |\n")
else:
    out.append('(pending)\n')
# 2. relaxation
rel = sorted([r for r in recs if r.get('batch') == 'relax'], key=lambda r: (r['d'], r['pD']))
out.append('\n## 2. Jacobian-quotient ("Hessian") relaxation ranks (exact mod p)\n')
out.append('Necessary condition B_D(N0) in J(F)_{2d} on the full E1-admissible family (J-type generators m F_i with lifts A m e_i, plus kernel reps with their lifts). rank Phi = #pairs means no nonzero solution over F_p-bar, hence over C.\n\n')
if rel:
    out.append('| d | parity | family (J + kernel) | #pairs | dim (S/J)_{2d} | rank Phi | excluded? |\n|---|---|---|---|---|---|---|\n')
    seen = set()
    for r in rel:
        k = (r['d'], r['pD'])
        if k in seen: continue
        seen.add(k)
        out.append(f"| {r['d']} | {'even' if r.get('pD',0)==0 else 'odd'} | {r.get('family_J','?')} + {r.get('family_kernel','?')} = {r.get('K')} | {r.get('npairs')} | {r.get('ring_dim_2d')} | {r.get('rank_Phi')} | {'**YES (proved, mod p -> C)**' if r.get('excluded_over_Fp') else 'no (rank-deficient by count)'} |\n")
else:
    out.append('(pending)\n')
# 3. searches
sr = [r for r in recs if r.get('batch') == 'e2search' and not str(r.get('tag','')).startswith('smoke')]
out.append('\n## 3. Full-E2 variable-projection multistart searches (normalized system, float, L-BFGS)\n')
if sr:
    out.append('| tag | d | parity | family dim K | n_c | rows_red | restarts | best rel | median rel | #rel<1e-3 | #rel<1e-6 | #rel<1e-8 | generic rank/nullity of T_D | time s |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n')
    agg = {}
    for r in sorted(sr, key=lambda r: (r['d'], r['pD'], r['tag'])):
        rs = r.get('rank_stats', {})
        out.append(f"| {r['tag']} | {r['d']} | {'even' if r['pD']==0 else 'odd'} | {r['K']} | {r['n_c']} | {r['rows_red']} | {r['restarts']} | {r['best_rel']:.3e} | {r['median_rel']:.3e} | {r['n_rel_lt_1e-3']} | {r['n_rel_lt_1e-6']} | {r['n_rel_lt_1e-8']} | {rs.get('rank')}/{rs.get('nullity')} | {r['time_s']:.0f} |\n")
        a = agg.setdefault((r['d'], r['pD']), {'restarts': 0, 'best': 1e9, 'K': r['K']})
        a['restarts'] += r['restarts']; a['best'] = min(a['best'], r['best_rel'])
    out.append('\nTotals per family: ' + '; '.join(f"d={k[0]} {'even' if k[1]==0 else 'odd'} (K={v['K']}): {v['restarts']} restarts, best rel {v['best']:.3e}" for k, v in sorted(agg.items())) + '\n')
else:
    out.append('(pending)\n')
# residual landscape from restarts files
out.append('\n### Residual landscape (per-restart files)\n')
for path in sorted(glob.glob(f'{MIR}/*.restarts.jsonl')):
    tag = os.path.basename(path).replace('.restarts.jsonl', '')
    if tag.startswith('smoke') or tag.startswith('planted'): continue
    rels = []
    for l in open(path):
        try: rels.append(json.loads(l)['rel'])
        except Exception: pass
    if not rels: continue
    rels.sort()
    q = lambda f: rels[min(len(rels)-1, int(f*len(rels)))]
    out.append(f"* `{tag}`: n={len(rels)}, min {rels[0]:.3e}, 10% {q(0.1):.3e}, median {q(0.5):.3e}, 90% {q(0.9):.3e}, max {rels[-1]:.3e}\n")
# 4. structured
stc = [r for r in recs if r.get('batch') == 'structured']
out.append('\n## 4. Structured families vs the E1-admissible space (exact mod p)\n')
if stc:
    out.append('dim(U cap W_d) with W_d = {D : A D in J} = J_d + kernel; "new" = part not already in J_d.\n\n| d | parity | family U | dim U | dim U cap W | dim U cap J | new |\n|---|---|---|---|---|---|---|\n')
    for r in sorted(stc, key=lambda r: (r['d'], r['parity'], -r['dim_U_cap_W'], r['family'])):
        if r['dim_U_cap_W'] == 0 and r['family'].startswith('vars[') and len(r['family']) < 12:
            continue  # 3-variable subsets with zero intersection are summarized
        out.append(f"| {r['d']} | {'even' if r['parity']==0 else 'odd'} | {r['family']} | {r['dim_U']} | {r['dim_U_cap_W']} | {r.get('dim_U_cap_J','')} | {r.get('dim_new','')} |\n")
    z3 = sum(1 for r in stc if r['dim_U_cap_W'] == 0 and r['family'].startswith('vars[') and len(r['family']) < 12)
    out.append(f'\n({z3} three-variable cone families with zero intersection omitted from the table.)\n')
else:
    out.append('(pending)\n')
ex = [r for r in recs if r.get('batch') == 'e2exact']
if ex:
    out.append('\n### Exact normalized-E2 tests of single denominators (mod p)\n\n| name | d | parity | e1 | consistent | rank | nullity | rows x cols | time s |\n|---|---|---|---|---|---|---|---|---|\n')
    for r in ex:
        out.append(f"| {r['name']} | {r['d']} | {r['pD']} | {r['e1'].get('e1') if isinstance(r['e1'], dict) else r['e1']} | {r.get('consistent')} | {r.get('rank')} | {r.get('nullity')} | {r.get('rows')}x{r.get('cols')} | {r.get('time_s')} |\n")
# 5. candidates
out.append('\n## 5. Candidates with residual < 1e-8 and exact confirmation attempts\n')
conf = [r for r in recs if r.get('batch') == 'confirm' or 'max_den' in r]
cands = []
for path in sorted(glob.glob(f'{MIR}/*.confirm.jsonl')):
    if 'planted' in path: continue
    for l in open(path):
        try: cands.append(json.loads(l))
        except Exception: pass
if cands:
    for c in cands:
        out.append(f"* {c.get('tag')}: max_den {c.get('max_den')} err {c.get('err'):.1e} -> {c.get('status')} (D: `{str(c.get('D'))[:150]}`)\n")
else:
    out.append('None so far: no restart reached rel < 1e-8 on this instance (see landscape above).\n')
# validation lines
val = [r for r in recs if r.get('batch') == 'validation']
out.append('\n## 6. Validation runs\n')
for r in val:
    out.append(f"* {r.get('family')}: {r.get('note')}\n")
out.append('\n## 7. Jobs\n')
for path in sorted(glob.glob(f'{MIR}/chain*.log')):
    out.append(f"* `{os.path.basename(path)}`: " + ' / '.join(l.strip() for l in open(path) if l.strip())[-600:] + '\n')
open(REP, 'w').write(''.join(out))
print('report rebuilt:', REP, len(recs), 'log lines')
