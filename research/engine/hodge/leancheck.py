"""Local (non-sandboxed) replica of the verifier compile sequence.

    python3 -m hodge.leancheck --instance block-q2 --dir /path/with/Witness.lean,Solution.lean,proof.json,witness.json [--jobs 2]

Builds HodgeTarget.lean / HodgeWitness.lean exactly as scripts/proposals/check.py does, then compiles
HodgeTarget, HodgeWitness, Proposal/<helpers>, Proposal/Solution, HodgeChecked with `lean -o` in a temp
layout, reporting wall time and peak memory for each step.
"""
from __future__ import annotations
import argparse, json, os, resource, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path('/home/kejace/devel/sextics')
sys.path.insert(0, str(ROOT / 'scripts/proposals'))
import yaml  # noqa
from check import target_source, witness_source, module_order  # noqa


def run(instance_id, folder, jobs=2, memory_mb=14000, keep=False):
    folder = Path(folder)
    files = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
    manifest = json.loads(files['proof.json'])
    instance = yaml.safe_load((ROOT / 'data' / instance_id / 'instance.rec').read_text())
    order, imports = module_order(files)
    lean_root = ROOT / 'formal/lean'
    lean_path = subprocess.check_output(['lake', 'env', 'printenv', 'LEAN_PATH'], cwd=lean_root).decode().strip()
    prefix = Path(subprocess.check_output(['lake', 'env', 'lean', '--print-prefix'], cwd=lean_root).decode().strip())
    lean = prefix / 'bin/lean'
    tmp = Path(tempfile.mkdtemp(prefix='hodge-leancheck-'))
    sources = tmp / 'source'; artifacts = tmp / 'artifacts'; protected = tmp / 'protected'; checked = tmp / 'checked'
    (sources / 'Proposal').mkdir(parents=True); (artifacts / 'Proposal').mkdir(parents=True); protected.mkdir(); checked.mkdir()
    for name, raw in files.items():
        if name.endswith('.lean'):
            (sources / 'Proposal' / name).write_bytes(raw)
    (sources / 'HodgeTarget.lean').write_text(target_source(instance))
    has_witness = 'witness.json' in files
    if has_witness:
        (sources / 'HodgeWitness.lean').write_text(witness_source(files['witness.json']))
    (sources / 'HodgeChecked.lean').write_text(
        'import ' + ('HodgeWitness' if has_witness else 'HodgeTarget') + '\nimport Proposal.Solution\n'
        'theorem HodgeSubmission.checked : HodgeComputer.HasCertificate '
        'HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q := ' + manifest['theorem'] + '\n' +
        ('theorem HodgeSubmission.checkedWitness : HodgeComputer.Verify '
         'HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q HodgeWitness.d '
         'HodgeWitness.D HodgeWitness.N HodgeWitness.Q := ' + manifest['witnessTheorem'] + '\n' if has_witness else ''))
    env = dict(os.environ)
    env['LEAN_PATH'] = lean_path + ':' + ':'.join(map(str, (protected, checked, artifacts)))
    protected_names = ['HodgeTarget.lean'] + (['HodgeWitness.lean'] if has_witness else [])
    sequence = [*protected_names, *order, 'HodgeChecked.lean']
    results = []
    total = 0.0
    for file in sequence:
        dest = protected if file in protected_names else checked if file == 'HodgeChecked.lean' else artifacts
        relative = 'Proposal/' + file if dest == artifacts else file
        out = dest / relative.replace('.lean', '.olean')
        t0 = time.time()
        r0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        p = subprocess.run([str(lean), '-j', str(jobs), '-M', str(memory_mb), '-o', str(out), relative],
                           cwd=sources, env=env, capture_output=True, text=True)
        dt = time.time() - t0
        total += dt
        r1 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        ok = p.returncode == 0
        results.append({'file': file, 'ok': ok, 'seconds': round(dt, 1), 'maxrss_mb': round(r1 / 1024)})
        print(f'[{file}] {"OK" if ok else "FAIL"} {dt:.1f}s', flush=True)
        if not ok:
            print(p.stdout[-6000:], p.stderr[-3000:], flush=True)
            break
    print(f'total {total:.1f}s; tmp dir {tmp}')
    if not keep and all(r['ok'] for r in results):
        import shutil; shutil.rmtree(tmp, ignore_errors=True)
    return results


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--instance', required=True)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--jobs', type=int, default=2)
    ap.add_argument('--keep', action='store_true')
    a = ap.parse_args()
    run(a.instance, a.dir, a.jobs, keep=a.keep)
