"""Post-processing: (1) verdicts for msolve leading-ideal outputs of the full relaxation systems (job B);
(2) summary tables of the wave-2 exclusion log; (3) writes packet README fragments.  Light; runs locally."""
import sys, os, json, glob, collections
sys.path.insert(0, '/tmp/hodge/results/klein-coupled-character/wave2-excl')
from w2lib import parse_msolve_leading, leading_ideal_irrelevant, LOG, PACKET

def msolve_verdicts(dirs):
    out = {}
    for d_ in dirs:
        for f in sorted(glob.glob(d_ + '/*.lead')):
            try:
                mons = parse_msolve_leading(f)
                ms = open(f[:-5] + '.ms').readline().strip().split(',')
                irr, info = leading_ideal_irrelevant(mons, len(ms))
                out[os.path.basename(f)] = {'n_vars': len(ms), 'n_lead': len(mons), 'irrelevant': irr,
                                            'missing': info['missing'][:10] if isinstance(info, dict) and not irr else None}
            except Exception as ex:
                out[os.path.basename(f)] = {'error': repr(ex)}
    return out

def summarize_log(path=LOG):
    rows = [json.loads(l) for l in open(path)] if os.path.exists(path) else []
    by = collections.defaultdict(list)
    for r in rows: by[r.get('exp')].append(r)
    lines = [f'total wave-2 log entries: {len(rows)}']
    for exp, rs in by.items():
        res = collections.Counter(str(r.get('result'))[:60] for r in rs)
        lines.append(f'- {exp}: n={len(rs)}; results: ' + '; '.join(f'{k} x{v}' for k, v in res.most_common(8)))
    return '\n'.join(lines), by

if __name__ == '__main__':
    dirs = sys.argv[1:] or glob.glob(PACKET + '/jobA_msolve') + glob.glob('/tmp/hodge/remote/klein-coupled-character/wave2-excl/jobA_msolve')
    v = msolve_verdicts(dirs)
    print(json.dumps(v, indent=1))
    s, by = summarize_log()
    print(s)
