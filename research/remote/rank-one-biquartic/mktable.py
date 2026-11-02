"""Render the experiment tables of the attack note from log.jsonl."""
import json, collections
rows = [json.loads(l) for l in open('/tmp/hodge/results/rank-one-biquartic/log.jsonl')]
fixed = [r for r in rows if r.get('family') != 'e1-kernel' and r.get('family') != 'timage']
kern = [r for r in rows if r.get('family') == 'e1-kernel']
tim = [r for r in rows if r.get('family') == 'timage']
out = []
out.append('### 2.1 Fixed-denominator tests (E1 = A·D ∈ J(F) mod p; then the full cleared linear system in (N,Q), compressed-rows rref mod p = 2^31−1)')
out.append('')
out.append('| D | deg | family | E1 | consistent | rank | nullity | rows×cols | s | note |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
seen = set()
for r in fixed:
    key = (r['D'], r.get('family'))
    if key in seen: continue
    seen.add(key)
    rc = f"{r.get('rows','')}×{r.get('cols','')}" if r.get('rows') else ''
    out.append(f"| `{r['D'][:60]}` | {r.get('d')} | {r.get('family')} | {r.get('e1')} | {r.get('consistent')} | {r.get('rank')} | {r.get('nullity')} | {rc} | {r.get('seconds')} | {(r.get('note') or '')[:80]} |")
out.append('')
out.append('### 2.2 E1-kernel (annihilator) sweeps: dim {D : T·D ∈ J}/J_d  (mod p)')
out.append('')
out.append('| target T | d | parity(D) | dim E1-space | dim J_d | kernel mod J | s |')
out.append('|---|---|---|---|---|---|---|')
seenk = set()
for r in kern:
    info = json.loads(r['note'])
    T = r['D'].split('target=')[1] if 'target=' in r['D'] else 'A = X0^3 Y0^3'
    key = (T, info['d'], info['parD'])
    if key in seenk: continue
    seenk.add(key)
    out.append(f"| {T} | {info['d']} | {info['parD']} | {info['dim_E1_space']} | {info['dim_J_d']} | **{info['kernel_mod_J']}** | {r.get('seconds')} |")
out.append('')
out.append('### 2.3 target_image: dimension of the space of Y-odd degree-6 targets served by D (18 = J_6 is the Lefschetz-trivial baseline)')
out.append('')
out.append('| D | deg | served dim | geometric (served − 18) | A served? | s |')
out.append('|---|---|---|---|---|---|')
for r in tim:
    out.append(f"| `{r['D']}` | {r.get('d')} | {r.get('rank')} | {r.get('nullity')} | {r.get('consistent')} | {r.get('seconds')} |")
print('\n'.join(out))
