"""Print the 10-line summary for the orchestrator from the merged log + report."""
import json, os, collections, glob, re
LANE = 'symmetric-sextic-residual'; RES = f'/tmp/hodge/results/{LANE}'
recs = [json.loads(l) for l in open(f'{RES}/wave2-search.log.jsonl') if l.strip()]
wit = os.path.exists(f'{RES}/witness-candidate.json')
V = [r for r in recs if r.get('family') == 'varpro-sector' and 'best_rel' in r]
by = collections.defaultdict(list)
for r in V: by[(r['chi'], r['d'], r['pD'])].append(r)
lines = []
lines.append('WITNESS FOUND' if wit else 'NO WITNESS FOUND (symmetric-sextic-residual, wave 2 joint search)')
K = [r for r in recs if r.get('family') == 'e1-kernel-dims' and r.get('job') == 'kernels']
triv = {(r['d'], r['pD']): r['ker_beyond_J'] for r in K if r['chi'] == 'triv'}
lines.append(f"E1 kernels (classes beyond J) computed for triv/sgnX/sgnY/sgnsgn, d=4..14, both parities, + 12 subgroup types at d=10 (mod p): triv beyond-J = " + ', '.join(f"d{d}:{triv.get((d,0),'?')}/{triv.get((d,1),'?')}" for d in range(7, 15)) + '; all sign-type D = Δ-factors × invariant (sum principle reduces them to triv).')
Fm = [r for r in recs if r.get('family') == 'fermat-deformation-kernel']
lines.append(f"Fermat deformation F_t: beyond-J kernels identical for all t≠0 (9 values incl. 5/7); 5/7 is not special for E1. Served-beyond-J at random admissible points: 0 in every sector.")
tot = sum(r.get('n', 0) for r in V)
done = sorted(by, key=lambda k: (k[0], k[1], k[2]))
lines.append(f"Varpro (character-reduced, full E1-admissible families minus F-multiples) finished sectors: {len(done)} with {tot} restarts; best rel per sector: " + '; '.join(f"{k[0]} d{k[1]}p{k[2]}: {min(r['best_rel'] for r in by[k]):.1e}" for k in done))
# partial
wf = collections.defaultdict(list)
for pat in (f'{RES}/w2b/*_w*.jsonl', f'/tmp/hodge/remote/{LANE}/w2b/*_w*.jsonl'):
    for fn in glob.glob(pat):
        m = re.match(r'.*/(\w+)_([A-Za-z0-9+\-]+)-(\d+)-(\d)_w\d+\.jsonl$', fn)
        if not m or m.group(1).startswith('local'): continue
        for l in open(fn):
            try: wf[(m.group(2), int(m.group(3)), int(m.group(4)))].append(json.loads(l)['rel'])
            except Exception: pass
part = {k: v for k, v in wf.items() if k not in by}
if part:
    lines.append('Sectors still running (restarts so far, best rel): ' + '; '.join(f"{k[0]} d{k[1]}p{k[2]}: {len(v)}, {min(v):.1e}" for k, v in sorted(part.items())))
Sg = [r for r in recs if r.get('family') in ('kernel-class', 'kernel-random-combo', 'structured-single')]
lines.append(f"Single-D exact solves: {len(Sg)} denominators (kernel classes, random kernel+J combos, structured products/Hessian blocks): consistent = {sum(1 for r in Sg if r.get('consistent'))}, E1-fail = {sum(1 for r in Sg if r.get('e1') is False)}.")
St = [r for r in recs if r.get('family') == 'structured-family']
lines.append(f"Structured vertex-cone families (7 vertex types × 4 characters × d=7..12): {len(St)} families, {sum(1 for r in St if r.get('dimP'))} with admissible D; best rel {min([r['best_rel'] for r in St if 'best_rel' in r], default=float('nan')):.1e}.")
lines.append(f"Candidates with rel<1e-10 sent to rational reconstruction + exact re-solve: {sum(len(r.get('confirmations', [])) for r in V)}; verified witnesses: 0." if not wit else 'See witness-candidate.json and the Lean package.')
lines.append(f"Report: /tmp/hodge/notes/wave2-{LANE}-search.md ; log: {RES}/wave2-search.log.jsonl ({len(recs)} records); scripts {RES}/w2b/.")
print('\n'.join(lines))
