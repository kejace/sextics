"""Fill placeholders from the final logs, build referee-report.md, copy the referee material into the packet (additive only), write the post files.
Usage: python3 assemble.py <SHA or 'SHA'>"""
import json, sys, os, shutil, collections
SHA = sys.argv[1] if len(sys.argv) > 1 else 'SHA'
W = '/tmp/hodge/work/referee42'; PK = '/tmp/hodge/results/four-plus-two-residual/wave2-excl'; POSTS = '/tmp/hodge/posts/four-plus-two-residual'
os.makedirs(POSTS, exist_ok=True)
def jl(path):
    out = []
    if not os.path.exists(path): return out
    for line in open(path):
        try: out.append(json.loads(line))
        except Exception: pass
    return out
sweep = [d for d in jl(f'{PK}/remote/excl-symsweep-local-k6-8.jsonl') if d.get('exp') == 'sym-fixedD']
cnt = collections.Counter((d['family'], d['result']) for d in sweep)
byk = collections.defaultdict(lambda: collections.Counter())
for (fam, res), n in cnt.items():
    k = int(fam.split('-')[0][3:]); kind = 'rep' if '-rep-' in fam else 'random'; par = fam[-1]
    byk[k][(kind, par, res)] += n
plan = {6: (20, 13), 7: (30, 30), 8: (41, 34)}
parts = []; allok = True
for k in sorted(byk):
    c = byk[k]
    reps0 = sum(n for (kind, par, res), n in c.items() if kind == 'rep' and par == '0'); reps1 = sum(n for (kind, par, res), n in c.items() if kind == 'rep' and par == '1')
    rnd = sum(n for (kind, par, res), n in c.items() if kind == 'random')
    bad = [(key, n) for key, n in c.items() if 'inconsistent over C [RANK]' not in key[2]]
    allok &= not bad
    parts.append(f"degree {k}: {reps0}/{plan[k][0]} even and {reps1}/{plan[k][1]} odd representatives plus {rnd} random same-parity combinations" + (f" (NON-RANK: {bad})" if bad else ''))
sweep_sentence = ('all of the following were inconsistent over $\\mathbb{C}$ by the nullity argument: ' if allok else 'results: ') + '; '.join(parts) + '.'
sweep_counts_report = '; '.join(parts) + ('; every record is "inconsistent over C [RANK]"' if allok else '; SOME RECORDS ARE NOT [RANK]')
t4 = jl(f'{W}/out/ref_t4.jsonl'); t5 = jl(f'{W}/out/ref_t5.jsonl')
def row(d):
    return f"{d.get('D')}: rank_p M = {d.get('rank_M_p')}, rank_p [M|b] = {d.get('rank_Mb_p')}, nullity_p = {d.get('nullity_p')}, exact nullity = {d.get('nullity_exact_Q')} -> {d.get('verdict')}"
fam_rows = [row(d) for d in t4 if d.get('exp') == 'fixedD-independent' and ('10-param' in str(d.get('D')) or '8-param' in str(d.get('D')))]
ker6_rows = [row(d) for d in t4 if d.get('exp') == 'fixedD-independent' and 'ker6' in str(d.get('D'))]
ker6_chk = [f"par{d['parity']}: {d['n_reps']} reps, rep0 in kernel exact = {d['A*rep0_in_J(H)_12_exact']}, not in J(H)_6 = {d['rep0_not_in_J(H)_6_exact']}" for d in t4 if d.get('exp') == 'ker6-rep-check']
t5_rows = [f"{d.get('D')} (sector-restricted, rows {d.get('rows')}, cols {d.get('cols')}): rank_p M = {d.get('rank_M_p')}, rank_p [M|b] = {d.get('rank_Mb_p')}, nullity_p = {d.get('nullity_p')}, exact nullity = {d.get('nullity_exact_Q')}, A*D in J(H) exact = {d.get('A*D_in_J(H)_exact')} -> {d.get('verdict')}" for d in t5 if d.get('exp') == 'symsweep-replay']
sweep_replay = ' | '.join(ker6_chk + ker6_rows + t5_rows) if (ker6_rows or t5_rows) else 'NO REPLAY COMPLETED'
family_results = ' | '.join(fam_rows) if fam_rows else 'NOT COMPLETED'
p1 = open(f'{W}/report-part1.md').read(); p2 = open(f'{W}/report-part2.md').read()
p2 = p2.replace('SWEEPREPLAY', sweep_replay).replace('SWEEPCOUNTS', sweep_counts_report).replace('FAMILYRESULTS', family_results)
report_pk = p1 + p2
report = p1 + p2 + f"\n## Snapshot\nPacket and referee material as of research-branch commit {SHA}.\n"
os.makedirs(f'{PK}/referee', exist_ok=True)
for fn in ('reflib.py', 'ref_t1.py', 'ref_t2.py', 'ref_t3.py', 'ref_t4.py', 'ref_t5.py', 'replay_d6jac2.py', 'assemble.py'):
    shutil.copy(f'{W}/{fn}', f'{PK}/referee/{fn}')
shutil.copytree(f'{W}/out', f'{PK}/referee/out', dirs_exist_ok=True)
open(f'{PK}/referee/referee-report.md', 'w').write(report_pk)
open(f'{POSTS}/referee-report.md', 'w').write(report)
c = open(f'{W}/comment-draft.md').read().replace('SWEEPCOUNTS', 'and the core-only kernel representatives of the exact E1 kernels tested one by one in the character-sector-restricted system, ' + sweep_sentence).replace('SHA', SHA)
open(f'{POSTS}/comment.md', 'w').write(c)
print('sweep:', sweep_sentence); print('replay:', sweep_replay[:500]); print('families:', family_results[:400]); print('report lines', report.count('\n'), '| comment lines', c.count('\n'), 'words', len(c.split()))
