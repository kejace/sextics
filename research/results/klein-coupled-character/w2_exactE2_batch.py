"""Run w2_exactE2.py for several (d:parD:k:prime) specs in parallel.  python3 -u w2_exactE2_batch.py procs timeout spec..."""
import sys, subprocess, time
from multiprocessing import Pool
procs = int(sys.argv[1]); timeout = int(sys.argv[2]); specs = sys.argv[3:]
def run(spec):
    d, parD, k, p = spec.split(':')
    t = time.time()
    r = subprocess.run(['python3', '-u', 'w2_exactE2.py', d, parD, k, 'J', str(timeout), p], capture_output=True, text=True,
                       cwd='/tmp/hodge/results/klein-coupled-character', env={**__import__('os').environ, 'OMP_NUM_THREADS': '1'})
    lines = [l for l in r.stdout.splitlines() if l.startswith(('system', 'singular'))]
    return spec, lines, round(time.time() - t, 1), r.stderr[-300:] if r.returncode else ''
if __name__ == '__main__':
    with Pool(procs) as pool:
        for spec, lines, s, err in pool.imap_unordered(run, specs):
            print(spec, lines, f'[{s}s]', err, flush=True)
    print('EXACTE2_BATCH_DONE', flush=True)
