"""sweep.py end-to-end on accepted witnesses: finds and writes exact witnesses, flags
perturbed/parity-mixed/over-cap candidates, resumes, and honors the timeout.
Scratch files are left under /tmp/hodge/results/_tests/."""
from common import *
import json, os, time
from hodge import core, sweep

def test_sweep_accepted():
    tmp = scratch_dir('sweep-')
    for iid in ['collision-h3', 'fermat-pair-rational-d4-strict']:
        F, A, q, _ = core.load_instance(iid)
        d, D, N, Q = core.load_witness(iid)
        cands = os.path.join(tmp, f'{iid}.jsonl')
        with open(cands, 'w') as f:
            f.write(json.dumps({'D': core.to_str(D), 'family': 'accepted'}) + '\n')
            D2 = dict(D); e = next(iter(D2)); D2[e] = D2[e] + 1
            f.write(json.dumps({'D': core.to_str(D2), 'family': 'perturbed'}) + '\n')
            f.write(core.to_str(D) + ' + X0*Y0' + ('*X1**%d' % (d - 2) if d > 2 else '') + '\n')  # mixed parity
            f.write(json.dumps({'D': f'({core.to_str(D)})**2', 'family': 'over-cap'}) + '\n')
        out = os.path.join(tmp, f'{iid}-out.jsonl'); log = os.path.join(tmp, f'{iid}-log.jsonl')
        wdir = os.path.join(tmp, f'{iid}-wit')
        t0 = time.time()
        res = sweep.sweep(iid, sweep.read_candidates(cands), out, workers=2, dmax=d, log_path=log,
                          witness_dir=wdir, timeout=600, verbose=False)
        assert len(res) == 4
        byfam = {r['family']: r for r in res}
        assert byfam['accepted'].get('witness') and byfam['accepted']['exact']
        assert byfam['perturbed']['consistent'] is False
        assert 'parity' in byfam['unspecified']['note']
        assert 'cap' in byfam['over-cap']['note']
        wfiles = os.listdir(wdir)
        assert len(wfiles) == 1
        w = json.load(open(os.path.join(wdir, wfiles[0])))
        Nw = [core.parse(s) for s in w['N']]
        ok, _ = core.verify_witness(F, A, q, core.parse(w['D']), Nw, core.parse(w['Q']))
        assert ok
        rec = json.loads(open(log).readline())
        for k in ('D', 'd', 'family', 'e1', 'consistent', 'rank', 'nullity', 'seconds', 'note'):
            assert k in rec
        # resume: nothing left to do
        res2 = sweep.sweep(iid, sweep.read_candidates(cands), out, workers=2, dmax=d, verbose=False)
        assert res2 == []
        print(f'  {iid}: sweep found/verified witness, controls flagged, resume ok ({time.time()-t0:.1f}s)')

def test_sweep_timeout_and_cli():
    import subprocess, sys
    tmp = scratch_dir('sweep-cli-')
    iid = 'hesse-triangle'
    d, D, N, Q = core.load_witness(iid)
    cands = os.path.join(tmp, 'c.jsonl')
    with open(cands, 'w') as f:
        f.write(json.dumps({'D': core.to_str(D), 'family': 'accepted'}) + '\n')
    out = os.path.join(tmp, 'out.jsonl')
    res = sweep.sweep(iid, sweep.read_candidates(cands), out, workers=1, timeout=1.0, verbose=False)
    assert len(res) == 1 and 'timeout' in res[0]['note'], res
    iid = 'collision-h3'
    d, D, N, Q = core.load_witness(iid)
    cands2 = os.path.join(tmp, 'c2.jsonl')
    with open(cands2, 'w') as f:
        f.write(core.to_str(D) + '\n')
    out2 = os.path.join(tmp, 'out2.jsonl')
    cp = subprocess.run([sys.executable, '-m', 'hodge.sweep', '--instance', iid, '--candidates', cands2,
                         '--out', out2, '--workers', '1', '--witness-dir', os.path.join(tmp, 'w'), '--symmetric'],
                        cwd=ROOT, capture_output=True, text=True, timeout=600)
    assert cp.returncode == 0, cp.stderr[-2000:]
    rec = json.loads(open(out2).readline())
    assert rec.get('witness') and rec['solver'] == 'symmetric'
    print('  timeout kill + CLI (--symmetric) ok')

if __name__ == '__main__':
    test_sweep_accepted(); test_sweep_timeout_and_cli()
    print('test_sweep OK')
