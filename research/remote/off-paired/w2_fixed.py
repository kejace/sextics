"""Step 4c: fixed structured candidates (semi-invariant, pure parity): sector E1 membership (exact basis, mod p)
then the complete sector-restricted fixed-D solve (compressed) with a wall-clock cap.
usage: python3 -u w2_fixed.py --cands <file.jsonl> --name <job> [--timeout 3600]
cands: {"D": "<expr>", "family": "..."} per line."""
import sys, os, json, time, argparse, signal
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, '/tmp/hodge/results/off-paired')
from fractions import Fraction
import w2lib as W
from hodge import core, symmetry, jacobian
ap = argparse.ArgumentParser(); ap.add_argument('--cands', required=True); ap.add_argument('--name', required=True)
ap.add_argument('--timeout', type=int, default=3600); ap.add_argument('--e1only', type=int, default=0)
args = ap.parse_args()
F, A, q = W.load(); G = symmetry.diagonal_group(F)
out = open(f'/tmp/hodge/results/off-paired/wave2-{args.name}.jsonl', 'a')
n = e1p = cons = 0; t00 = time.time()


class TO(Exception):
    pass


def handler(signum, frame):
    raise TO()


for line in open(args.cands):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    cand = json.loads(line); t0 = time.time()
    D = core.parse(cand['D']); d = core.degree(D); c = W.poly_chi(D); par = core.parity_of(D)
    rec = {'D': cand['D'][:300], 'family': cand.get('family'), 'd': d, 'chi': list(c) if c else None, 'parity': par}
    n += 1
    if c is None or par is None:
        rec['note'] = 'not semi-invariant / mixed parity'; out.write(json.dumps(rec) + '\n'); out.flush(); continue
    divF, _ = jacobian.divisible_by(F, D)
    rec['F_divides'] = bool(divF)
    S = W.SectorSpaces(F, A, q, d, c)
    Dm = S.Dmons; K = len(Dm)
    e1rows = [S._modrow([Fraction(v) for v in row]) for row in S.E1]
    vrow = S._modrow(S._vec(D, Dm))
    e1 = S._rank_rows(e1rows + [vrow], K) == S._rank_rows(e1rows, K)
    rec['e1'] = bool(e1); rec['dimE1_sector'] = S.dimE1
    if e1:
        e1p += 1
        if not args.e1only and not divF:
            signal.signal(signal.SIGALRM, handler); signal.alarm(args.timeout)
            try:
                st, data = symmetry.solve_fixed_D_symmetric(F, A, q, D, group=G, compressed=True, max_primes=3)
                rec['status'] = st; rec.update({k: data.get(k) for k in ('rank', 'nullity', 'rows', 'cols')})
                if st == 'solved':
                    cons += 1
                    wit = {'d': d, 'D': cand['D'], 'N': [W.D_str(x) for x in data['N']], 'Q': W.D_str(data['Q'])}
                    json.dump(wit, open('/tmp/hodge/results/off-paired/witness-candidate.json', 'w'), indent=1)
                    W.log_batch({'step': 'WITNESS', 'name': args.name, 'd': d, 'D': cand['D'][:300]})
                    print('WITNESS FOUND', cand['D'][:200], flush=True)
            except TO:
                rec['status'] = 'timeout'
            except Exception as ex:
                rec['status'] = f'error: {ex}'
            finally:
                signal.alarm(0)
    rec['seconds'] = round(time.time() - t0, 1)
    out.write(json.dumps(rec) + '\n'); out.flush()
    print(json.dumps({k: rec[k] for k in rec if k != 'D'}), cand['D'][:80], flush=True)
W.log_batch({'step': 'sweep', 'name': args.name, 'family': os.path.basename(args.cands), 'n': n, 'e1_pass': e1p,
             'consistent': cons, 'seconds': round(time.time() - t00, 1), 'degrees': 'see wave2-%s.jsonl' % args.name})
os.system('python3 /tmp/hodge/results/off-paired/w2_report.py > /dev/null 2>&1')
